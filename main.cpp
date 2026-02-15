// ============================================================================
// main.cpp - Simulador N-body con integración híbrida (KS + campo)
// ============================================================================
// 
// EVOLUCIÓN DEL CÓDIGO:
// ---------------------
// Versión 1.0 (Original): Simulación básica de binaria con integrador híbrido
// Versión 1.1: Se añadió logger de regímenes (regimes.csv)
// Versión 1.2: Se corrigió la interfaz de Integrator (ahora usa vector<bool> used)
// Versión 1.3: Se implementó selección óptima de binarias (no greedy)
// Versión 1.4: Se añadió time-centering (drift-kick) en HybridIntegrator
// Versión 1.5: Se completó KSIntegrator con subciclado
// Versión 1.6: Se añadieron condiciones iniciales (InitialConditions)
// Versión 1.7: Se implementó línea de comandos para seleccionar escenario
// Versión 1.8: Se añadió documentación extensa y modo debug
//
// ============================================================================
// USO DESDE LÍNEA DE COMANDOS:
// ============================================================================
// nbody_sim.exe [escenario] [opciones]
//
// Escenarios disponibles:
//   binary   - Binaria Kepleriana circular (default)
//              Masas: 1.0, 1.0 | Separación: 1.0 | Excentricidad: 0.0
//              Propósito: Testear KS y conservación de energía
//
//   field    - Binaria + estrella de campo
//              Binaria: masas 1.0, separación 1.0
//              Campo: masa 1.0 a distancia 10
//              Propósito: Testear el integrador híbrido (KS + campo)
//
//   figure8  - Órbita figura-8 de 3 cuerpos iguales
//              Configuración de Chenciner & Montgomery (2000)
//              Propósito: Testear sistemas multi-cuerpo estables
//
//   solar    - Sistema solar simplificado (Sol + Tierra)
//              Sol: masa 1.0 en el centro
//              Tierra: masa 3e-6, órbita circular a 1 UA
//              Propósito: Escala realista de masas
//
//   random   - Sistema aleatorio de N cuerpos
//              Por defecto: 10 cuerpos en radio 5.0
//              Propósito: Tests estadísticos
//
//   plummer  - Cúmulo de Plummer (distribución realista)
//              Por defecto: 100 cuerpos
//              Propósito: Simulaciones astrofísicas realistas
//
// ============================================================================
// ARCHIVOS DE SALIDA:
// ============================================================================
// - simulation.csv : Datos de energía, momento, momento angular por paso
// - regimes.csv    : Registro de activaciones/desactivaciones de KS
//
// ============================================================================
// COMPILACIÓN (con CMake):
// ============================================================================
// mkdir build && cd build
// cmake .. -G "Ninja"
// cmake --build .
//
// ============================================================================

#include <iostream>
#include <memory>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>

// Includes del proyecto
#include "nbody_system.h"
#include "euler_integrator.h"
#include "leapfrog_integrator.h"
#include "rk4_integrator.h"
#include "velocity_verlet_integrator.h"
#include "rk45_integrator.h"
#include "hybrid_integrator.h"
#include "data_logger.h"
#include "regime_logger.h"
#include "initial_conditions.h"
#include "system_analyzer.h"

// ============================================================================
// PROTOTIPOS DE FUNCIONES AUXILIARES
// ============================================================================
void print_usage();
void print_system_info(const NBodySystem& system, const std::string& scenario);
void print_progress(int step, int total, double E, double dE);
void print_final_stats(const NBodySystem& system, double E0, int steps);

// ============================================================================
// FUNCIÓN PRINCIPAL
// ============================================================================
int main(int argc, char* argv[]) {
    // ------------------------------------------------------------------------
    // 1. CONFIGURACIÓN DESDE LÍNEA DE COMANDOS
    // ------------------------------------------------------------------------
    // Valores por defecto
    std::string scenario = "binary";
    double dt = 0.01;
    int steps = 1000;
    double r_close = 0.5;
    double ks_dt = 1e-4;
    bool verbose = false;
    
    // Procesar argumentos (sintaxis simple: escenario [dt] [steps])
    if (argc >= 2) {
        scenario = argv[1];
    }
    if (argc >= 3) {
        dt = std::stod(argv[2]);
    }
    if (argc >= 4) {
        steps = std::stoi(argv[3]);
    }
    if (argc >= 5) {
        r_close = std::stod(argv[4]);
    }
    if (argc >= 6) {
        ks_dt = std::stod(argv[5]);
    }
    
    // ------------------------------------------------------------------------
    // 2. SELECCIÓN DEL ESCENARIO (condiciones iniciales)
    // ------------------------------------------------------------------------
    NBodySystem system;
    
    try {
        if (scenario == "binary") {
            system = InitialConditions::kepler_binary(1.0, 0.0);
        }
        else if (scenario == "field") {
            system = InitialConditions::binary_with_field();
        }
        else if (scenario == "figure8") {
            system = InitialConditions::figure_eight();
        }
        else if (scenario == "solar") {
            system = InitialConditions::solar_system();
        }
        else if (scenario == "random") {
            system = InitialConditions::random_system(10, 5.0);
        }
        else if (scenario == "plummer") {
            system = InitialConditions::plummer_cluster(100);
        }
        else {
            std::cerr << "\n❌ Error: escenario '" << scenario << "' no reconocido\n\n";
            print_usage();
            return 1;
        }
    }
    catch (const std::exception& e) {
        std::cerr << "\n❌ Error al crear condiciones iniciales: " << e.what() << "\n";
        return 1;
    }
    
    // ------------------------------------------------------------------------
    // 3. MOSTRAR INFORMACIÓN DEL SISTEMA
    // ------------------------------------------------------------------------
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "🚀 SIMULADOR N-BODY CON INTEGRACIÓN HÍBRIDA\n";
    std::cout << std::string(60, '=') << "\n";
    std::cout << "Escenario        : " << scenario << "\n";
    std::cout << "Número de cuerpos: " << system.bodies.size() << "\n";
    std::cout << "Paso de tiempo   : " << dt << "\n";
    std::cout << "Número de pasos  : " << steps << "\n";
    std::cout << "Radio crítico    : " << r_close << "\n";
    std::cout << "Paso interno KS  : " << ks_dt << "\n";
    std::cout << std::string(60, '-') << "\n";
    
    // Mostrar información detallada de cada cuerpo (opcional)
    if (verbose) {
        print_system_info(system, scenario);
    }
    
    // ------------------------------------------------------------------------
    // 4. LOGGERS (archivos de salida)
    // ------------------------------------------------------------------------
    // Logger de regímenes (activación de KS)
    auto regime_logger = std::make_unique<RegimeLogger>("regimes.csv");
    std::cout << "📊 Archivo de regímenes: regimes.csv\n";
    
    // Logger de datos (energía, momento, etc.)
    DataLogger data_logger("simulation.csv");
    std::cout << "📈 Archivo de simulación: simulation.csv\n";
    std::cout << std::string(60, '-') << "\n";
    
    // ------------------------------------------------------------------------
    // 5. INVARIANTES INICIALES (para calcular errores)
    // ------------------------------------------------------------------------
    const double E0 = system.total_energy();
    const Vec3 P0 = system.total_momentum();
    const Vec3 L0 = system.total_angular_momentum();
    
    std::cout << "Energía inicial   : " << std::setprecision(12) << E0 << "\n";
    std::cout << "Momento inicial   : (" << P0.x << ", " << P0.y << ", " << P0.z << ")\n";
    std::cout << "Mom. angular init.: (" << L0.x << ", " << L0.y << ", " << L0.z << ")\n";
    std::cout << std::string(60, '=') << "\n\n";
    
    // ------------------------------------------------------------------------
    // 6. CONFIGURACIÓN DEL INTEGRADOR HÍBRIDO
    // ------------------------------------------------------------------------
    // El integrador híbrido combina:
    // - VelocityVerlet para el campo (fondo)
    // - KS para binarias cerradas
    // - Time-centering (drift-kick) para consistencia temporal
    // - Selección óptima de binarias (por energía)
    auto integrator = std::make_unique<HybridIntegrator>(
        std::make_unique<VelocityVerletIntegrator>(),  // integrador de fondo
        r_close,                                         // radio de detección
        ks_dt,                                           // paso interno KS
        regime_logger.get()                              // logger opcional
    );
    
    // ------------------------------------------------------------------------
    // 7. BUCLE PRINCIPAL DE SIMULACIÓN
    // ------------------------------------------------------------------------
    std::cout << "⏳ Simulando...\n";
    std::cout << "step      Energía         dE/E\n";
    std::cout << std::string(60, '-') << "\n";
    
    for (int step = 0; step < steps; ++step) {
        // Todos los cuerpos están activos (ninguno pre-marcado)
        std::vector<bool> used(system.bodies.size(), false);
        
        // Un paso del integrador híbrido
        integrator->step(system, dt, used);
        
        // Calcular invariantes actuales
        const double E = system.total_energy();
        const Vec3 P = system.total_momentum();
        const Vec3 L = system.total_angular_momentum();
        
        // Guardar en archivo CSV
        data_logger.log(
            step,
            E, E - E0,
            norm(P), norm(P) - norm(P0),
            norm(L), norm(L) - norm(L0)
        );
        
        // Mostrar progreso (cada 100 pasos o si hay warning)
        if (step % 100 == 0 || step == steps-1) {
            double dE_rel = std::abs(E - E0) / std::abs(E0);
            std::cout << std::setw(4) << step << "   "
                      << std::setprecision(10) << E << "   "
                      << std::setprecision(4) << dE_rel << "\n";
            
            // Warning si el error es grande
            if (dE_rel > 1e-6 && step > 0) {
                std::cout << "⚠️  Warning: Error de energía > 1e-6 en paso " << step << "\n";
            }
        }
    }
    
    // ------------------------------------------------------------------------
    // 8. ESTADÍSTICAS FINALES
    // ------------------------------------------------------------------------
    std::cout << std::string(60, '=') << "\n";
    std::cout << "✅ SIMULACIÓN COMPLETADA\n";
    std::cout << std::string(60, '-') << "\n";
    
    double E_final = system.total_energy();
    double dE_rel = std::abs(E_final - E0) / std::abs(E0);
    
    std::cout << "Energía final     : " << std::setprecision(12) << E_final << "\n";
    std::cout << "Error relativo    : " << std::setprecision(4) << dE_rel;
    
    if (dE_rel < 1e-10) std::cout << " (🔥 Excelente!)\n";
    else if (dE_rel < 1e-8) std::cout << " (👍 Muy bueno)\n";
    else if (dE_rel < 1e-6) std::cout << " (⚠️  Aceptable)\n";
    else std::cout << " (❌ Malo - revisa)\n";
    
    std::cout << std::string(60, '=') << "\n\n";
    
    // Análisis adicional con SystemAnalyzer
    SystemAnalyzer analyzer;
    auto snapshot = analyzer.analyze(system, steps * dt);
    
    std::cout << "📊 ANÁLISIS DEL SISTEMA FINAL:\n";
    std::cout << "  Binarias detectadas: " << snapshot.n_binaries << "\n";
    std::cout << "  Separación mínima  : " << snapshot.min_separation << "\n";
    std::cout << "  Separación máxima  : " << snapshot.max_separation << "\n";
    std::cout << "  Radio de media masa: " << analyzer.half_mass_radius(system) << "\n";
    std::cout << std::string(60, '=') << "\n";
    
    return 0;
}

// ============================================================================
// FUNCIONES AUXILIARES
// ============================================================================

void print_usage() {
    std::cout << "\n📋 USO: nbody_sim.exe [escenario] [dt] [steps] [r_close] [ks_dt]\n";
    std::cout << "\nEscenarios disponibles:\n";
    std::cout << "  binary   - Binaria Kepleriana circular (default)\n";
    std::cout << "  field    - Binaria + estrella de campo\n";
    std::cout << "  figure8  - Órbita figura-8 de 3 cuerpos\n";
    std::cout << "  solar    - Sistema solar simplificado\n";
    std::cout << "  random   - Sistema aleatorio (10 cuerpos)\n";
    std::cout << "  plummer  - Cúmulo de Plummer (100 cuerpos)\n";
    std::cout << "\nEjemplos:\n";
    std::cout << "  nbody_sim.exe binary\n";
    std::cout << "  nbody_sim.exe field 0.01 1000\n";
    std::cout << "  nbody_sim.exe figure8 0.005 5000 0.3 1e-5\n";
    std::cout << std::endl;
}

void print_system_info(const NBodySystem& system, const std::string& scenario) {
    std::cout << "\n📋 DETALLE DE CUERPOS:\n";
    for (size_t i = 0; i < system.bodies.size(); ++i) {
        const auto& b = system.bodies[i];
        std::cout << "  Cuerpo " << i << ": masa=" << b.mass
                  << "  pos=(" << b.position.x << ", " << b.position.y << ", " << b.position.z << ")"
                  << "  vel=(" << b.velocity.x << ", " << b.velocity.y << ", " << b.velocity.z << ")\n";
    }
    std::cout << std::endl;
}

void print_progress(int step, int total, double E, double dE) {
    // Función simplificada - ya se hace en el bucle principal
}

void print_final_stats(const NBodySystem& system, double E0, int steps) {
    // Función simplificada - ya se hace al final
}

/*
============================================================================
main.cpp - Simulador N-body con integración híbrida (KS + campo)
============================================================================

#include <iostream>
#include <memory>
#include <iomanip>

#include "nbody_system.h"
#include "euler_integrator.h"
#include "leapfrog_integrator.h"
#include "rk4_integrator.h"
#include "velocity_verlet_integrator.h"
#include "rk45_integrator.h"
#include "hybrid_integrator.h"
#include "data_logger.h"
#include "regime_logger.h"  // Añadido para el logger opcional

int main() {
    // =========================
    // Construcción del sistema
    // =========================
    NBodySystem system;

    system.bodies.push_back({{-1, 0, 0}, {0, -0.5, 0}, 1.0});
    system.bodies.push_back({{ 1, 0, 0}, {0,  0.5, 0}, 1.0});

    // =========================
    // Parámetros de simulación
    // =========================
    double dt = 0.01;
    const int steps = 1000;

    // =========================
    // Logger de regímenes (opcional)
    // =========================
    auto regime_logger = std::make_unique<RegimeLogger>("regimes.csv");

    // =========================
    // Selección del integrador
    // =========================
    std::unique_ptr<Integrator> integrator;

    // Versión actualizada con la nueva firma:
    // HybridIntegrator( far_integrator, r_close, ks_internal_dt, logger )
    integrator = std::make_unique<HybridIntegrator>(
        std::make_unique<VelocityVerletIntegrator>(),  // integrador de fondo
        0.5,                                            // r_close
        1e-4,                                           // ks_internal_dt (paso interno pequeño)
        regime_logger.get()                             // logger (opcional)
    );

    // =========================
    // Invariantes iniciales
    // =========================
    const double E0 = system.total_energy();
    const Vec3 P0 = system.total_momentum();
    const Vec3 L0 = system.total_angular_momentum();

    // =========================
    // Logger de datos
    // =========================
    DataLogger data_logger("simulation.csv");

    std::cout << std::setprecision(12);
    std::cout << "# step E dE |P| d|P| |L| d|L|\n";

    // =========================
    // Bucle principal
    // =========================
    for (int step = 0; step < steps; ++step) {
        // Vector used vacío porque no hay cuerpos pre-integrados
        std::vector<bool> used(system.bodies.size(), false);
        integrator->step(system, dt, used);

        const double E = system.total_energy();
        const Vec3 P = system.total_momentum();
        const Vec3 L = system.total_angular_momentum();

        data_logger.log(
            step,
            E, E - E0,
            norm(P), norm(P) - norm(P0),
            norm(L), norm(L) - norm(L0)
        );

        std::cout
            << step << " "
            << E << " " << (E - E0) << " "
            << norm(P) << " " << (norm(P) - norm(P0)) << " "
            << norm(L) << " " << (norm(L) - norm(L0))
            << "\n";
    }

    return 0;
}

*/