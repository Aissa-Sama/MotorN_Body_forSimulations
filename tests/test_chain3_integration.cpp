// tests/test_chain3_integration.cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include "nbody_system.h"
#include "chain3_integrator.h"
#include "initial_conditions.h"

/**
 * @brief Test de integración de la cadena para 3 cuerpos.
 * 
 * Verifica que la cadena pueda integrar el sistema figura-8
 * con la versión simplificada (sin control de error, sin tiempo físico).
 * 
 * ============================================================================
 * CONFIGURACIÓN ACTUAL (MODO DEPURACIÓN)
 * ============================================================================
 * - Control adaptativo: DESACTIVADO
 * - Tiempo físico: DESACTIVADO
 * - Número de pasos: 10000 fijos
 * - Paso ficticio (dtau): 1e-3
 * - Energía: NO se usa en las ecuaciones
 * ============================================================================
 */

// ----------------------------------------------------------------------------
// Función para calcular el período aproximado de la figura-8
// ----------------------------------------------------------------------------
double figure8_period() {
    return 6.3259;  // Período conocido de la órbita figura-8
}

// ----------------------------------------------------------------------------
// Función para calcular el error relativo
// ----------------------------------------------------------------------------
double rel_error(double computed, double reference) {
    return std::abs(computed - reference) / (std::abs(reference) + 1e-15);
}

// ----------------------------------------------------------------------------
// Función para calcular la diferencia angular entre vectores
// ----------------------------------------------------------------------------
double angular_difference(const Vec3& a, const Vec3& b) {
    double dot_product = a.dot(b);
    double norms = a.norm() * b.norm();
    if (norms < 1e-15) return 0.0;
    double cos_angle = dot_product / norms;
    // Asegurar que está en el rango [-1, 1] por errores numéricos
    if (cos_angle > 1.0) cos_angle = 1.0;
    if (cos_angle < -1.0) cos_angle = -1.0;
    return std::acos(cos_angle);
}

int main() {
    // ------------------------------------------------------------------------
    // 1. INICIALIZACIÓN Y CONFIGURACIÓN
    // ------------------------------------------------------------------------
    std::cout << std::setprecision(15);
    std::cout << "=== Test de integración de cadena (3 cuerpos) ===\n";
    std::cout << "=== Período de la figura-8 ≈ " << figure8_period() << " ===\n\n";
    
    // Crear sistema figura-8
    NBodySystem system = InitialConditions::figure_eight();
    double E0 = system.total_energy();
    Vec3 P0 = system.total_momentum();
    Vec3 L0 = system.total_angular_momentum();
    
    // Guardar posiciones iniciales para referencia
    std::vector<Vec3> r0;
    for (const auto& b : system.bodies) {
        r0.push_back(b.position);
    }
    
    std::cout << "Configuración:\n";
    std::cout << "  Tiempo objetivo: " << figure8_period() << "\n";
    std::cout << "  Paso ficticio interno: " << 1e-3 << "\n";
    std::cout << "  Control adaptativo: NO\n";
    std::cout << "  Tiempo físico: NO\n";
    std::cout << "  Número de pasos: 10000\n\n";
    
    std::cout << "Sistema inicial (figura-8):\n";
    for (size_t i = 0; i < system.bodies.size(); ++i) {
        const auto& b = system.bodies[i];
        std::cout << "  Cuerpo " << i << ": masa=" << b.mass
                  << "  pos=(" << b.position.x << ", " << b.position.y << ", " << b.position.z << ")"
                  << "  vel=(" << b.velocity.x << ", " << b.velocity.y << ", " << b.velocity.z << ")\n";
    }
    std::cout << "  Energía total inicial: " << E0 << "\n";
    std::cout << "  Momento angular total: (" << L0.x << ", " << L0.y << ", " << L0.z << ")\n\n";

    // ------------------------------------------------------------------------
    // 2. INICIALIZAR CADENA (transformar a variables KS)
    // ------------------------------------------------------------------------
    // PARÁMETRO: internal_dt - paso ficticio interno
    // En modo depuración usamos 1e-3 (1000 pasos por unidad de τ)
    Chain3Integrator chain3(1e-3);
    
    Chain3State state = chain3.initialize(system, 0, 1, 2);
    
    std::cout << "Estado de cadena inicial:\n";
    std::cout << "  |Q1|² = " << state.Q1.norm2() << " (debe ser ≈ |r12|)\n";
    std::cout << "  |Q2|² = " << state.Q2.norm2() << " (debe ser ≈ |r23|)\n";
    std::cout << "  Energía interna = " << state.energy << "\n";
    std::cout << "  CM: pos=(" << state.cm_pos.x << ", " << state.cm_pos.y << ", " << state.cm_pos.z << ")\n";
    std::cout << "  CM: vel=(" << state.cm_vel.x << ", " << state.cm_vel.y << ", " << state.cm_vel.z << ")\n\n";

    // ------------------------------------------------------------------------
    // 3. CONFIGURAR PARÁMETROS DE INTEGRACIÓN (no se usan en modo depuración)
    // ------------------------------------------------------------------------
    IntegrationParams params;
    params.rel_tol = 1e-8;      // No se usa
    params.abs_tol = 1e-12;     // No se usa
    params.min_dtau = 1e-12;    // No se usa
    params.max_dtau = 1.0;       // No se usa
    params.safety_factor = 0.9;  // No se usa
    params.verbose = true;       // Mostrar progreso

    // ------------------------------------------------------------------------
    // 4. EJECUTAR INTEGRACIÓN SIMPLIFICADA
    // ------------------------------------------------------------------------
    double t_target = figure8_period();  // No se usa realmente
    double t_achieved;
    
    std::cout << "Integrando " << 10000 << " pasos con dtau = " << 1e-3 << "...\n";
    std::cout << "----------------------------------------\n";
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    chain3.integrate(state, t_target, t_achieved, params, system);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration<double>(end_time - start_time).count();
    
    std::cout << "----------------------------------------\n";
    std::cout << "Tiempo alcanzado: " << t_achieved << " (no se usa en modo depuración)\n";
    std::cout << "Tiempo de ejecución: " << elapsed << " segundos\n\n";

    // ------------------------------------------------------------------------
    // 5. RECONSTRUIR SISTEMA FÍSICO
    // ------------------------------------------------------------------------
    NBodySystem system_final;
    system_final.bodies.resize(3);
    system_final.bodies[0].mass = state.m1();
    system_final.bodies[1].mass = state.m2();
    system_final.bodies[2].mass = state.m3();
    chain3.write_back(state, system_final, 0, 1, 2);
    
    std::cout << "Sistema final (reconstruido):\n";
    for (size_t i = 0; i < system_final.bodies.size(); ++i) {
        const auto& b = system_final.bodies[i];
        std::cout << "  Cuerpo " << i << ": pos=(" 
                  << b.position.x << ", " << b.position.y << ", " << b.position.z << ")"
                  << "  vel=(" << b.velocity.x << ", " << b.velocity.y << ", " << b.velocity.z << ")\n";
    }
    std::cout << "\n";

    // ------------------------------------------------------------------------
    // 6. VERIFICAR INVARIANTES
    // ------------------------------------------------------------------------
    double Ef = system_final.total_energy();
    Vec3 Pf = system_final.total_momentum();
    Vec3 Lf = system_final.total_angular_momentum();
    
    // Calcular errores relativos
    double dE_rel = rel_error(Ef, E0);
    double dP_rel = rel_error(Pf.norm(), P0.norm());
    double dL_rel = rel_error(Lf.norm(), L0.norm());
    
    // Calcular error promedio en posición (comparando distancias, no posiciones absolutas)
    double pos_error_sum = 0.0;
    for (int i = 0; i < 3; ++i) {
        // Comparar distancias al centro de masa, no posiciones absolutas
        Vec3 r_rel_orig = system.bodies[i].position - state.cm_pos;
        Vec3 r_rel_final = system_final.bodies[i].position - state.cm_pos;
        pos_error_sum += (r_rel_final - r_rel_orig).norm();
    }
    double pos_error_avg = pos_error_sum / 3.0;
    
    std::cout << "==================================================\n";
    std::cout << "RESULTADOS DE LA VERIFICACIÓN:\n";
    std::cout << "==================================================\n";
    std::cout << "Energía final: " << Ef << "\n";
    std::cout << "Error relativo en energía: " << dE_rel << "\n";
    std::cout << "Error relativo en momento: " << dP_rel << "\n";
    std::cout << "Error relativo en momento angular: " << dL_rel << "\n";
    std::cout << "Error promedio en posición: " << pos_error_avg << "\n";
    std::cout << "--------------------------------------------------\n";
    
    // ------------------------------------------------------------------------
    // 7. CRITERIOS DE ÉXITO (más flexibles para modo depuración)
    // ------------------------------------------------------------------------
    bool energia_ok = (dE_rel < 1e-2);      // 1% de error en energía
    bool momento_ok = (dP_rel < 1e-12);      // Momento debe conservarse exactamente
    bool angular_ok = (dL_rel < 1e-2);       // 1% de error en momento angular
    bool posicion_ok = (pos_error_avg < 1.0); // Error de posición < 1 unidad
    
    std::cout << "\nCRITERIOS DE ÉXITO (modo depuración):\n";
    std::cout << "  Conservación de energía ( < 1e-2 ): " 
              << (energia_ok ? "✅" : "❌") << " (" << dE_rel << ")\n";
    std::cout << "  Conservación de momento ( < 1e-12): " 
              << (momento_ok ? "✅" : "❌") << " (" << dP_rel << ")\n";
    std::cout << "  Conservación de mom. angular ( < 1e-2 ): " 
              << (angular_ok ? "✅" : "❌") << " (" << dL_rel << ")\n";
    std::cout << "  Precisión de posición ( < 1.0 ): " 
              << (posicion_ok ? "✅" : "❌") << " (" << pos_error_avg << ")\n";
    
    std::cout << "\n==================================================\n";
    
    int pruebas_pasadas = 0;
    if (energia_ok) pruebas_pasadas++;
    if (momento_ok) pruebas_pasadas++;
    if (angular_ok) pruebas_pasadas++;
    if (posicion_ok) pruebas_pasadas++;
    
    std::cout << "RESULTADO FINAL: " << pruebas_pasadas << "/4 pruebas pasadas\n";
    
    if (pruebas_pasadas >= 3) {
        std::cout << "✅ El integrador funciona aceptablemente en modo depuración\n";
        return 0;
    } else if (pruebas_pasadas >= 2) {
        std::cout << "⚠️ El integrador funciona parcialmente - revisar implementación\n";
        return 1;
    } else {
        std::cout << "❌ El integrador NO funciona - revisar implementación\n";
        return 2;
    }
}