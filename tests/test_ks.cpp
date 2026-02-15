// tests/test_ks.cpp
#include <iostream>
#include <cmath>
#include <iomanip>
#include "binary_state.h"
#include "ks_integrator.h"
#include "vec3.h"

// Constante PI portable
const double PI = 3.14159265358979323846;

// Función para calcular el período Kepleriano
double kepler_period(double a, double M_total) {
    return 2.0 * PI * std::sqrt(a * a * a / M_total);
}

bool test_ks_conservation() {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "🔬 TEST: Conservación de KS en órbita circular\n";
    std::cout << std::string(60, '=') << "\n";
    
    // Crear una binaria Kepleriana circular
    double m1 = 1.0, m2 = 1.0;
    double a = 1.0;  // semieje mayor
    double M = m1 + m2;
    double v = std::sqrt(M / a) * 0.5;  // velocidad para cada cuerpo
    
    Body body1{{-a/2, 0, 0}, {0, -v, 0}, m1};
    Body body2{{ a/2, 0, 0}, {0,  v, 0}, m2};
    
    BinaryState binary(body1, body2);
    
    // Calcular invariantes iniciales
    double r0 = binary.separation();
    Vec3 v0 = binary.relative_velocity();
    double mu = binary.reduced_mass();
    double E0 = 0.5 * mu * dot(v0, v0) - (m1 * m2) / r0;
    double T = kepler_period(a, M);
    
    std::cout << "Parámetros iniciales:\n";
    std::cout << "  Semieje mayor a = " << a << "\n";
    std::cout << "  Período T = " << T << "\n";
    std::cout << "  Energía E0 = " << std::setprecision(12) << E0 << "\n";
    std::cout << "  Separación r0 = " << r0 << "\n";
    
    // Probar diferentes pasos internos
    double ks_dt_values[] = {1e-3, 1e-4, 1e-5};
    
    for (double ks_dt : ks_dt_values) {
        std::cout << "\n" << std::string(40, '-') << "\n";
        std::cout << "Paso interno KS = " << ks_dt << "\n";
        
        KSIntegrator ks(ks_dt);
        BinaryState binary_copy = binary;  // Reiniciar cada vez
        
        // Integrar por un período
        ks.integrate(binary_copy, T);
        
        // Verificar después de un período
        double r1 = binary_copy.separation();
        Vec3 v1 = binary_copy.relative_velocity();
        double E1 = 0.5 * mu * dot(v1, v1) - (m1 * m2) / r1;
        
        double dE = std::abs(E1 - E0) / std::abs(E0);
        double dr = std::abs(r1 - r0) / r0;
        
        std::cout << "  Energía final: " << std::setprecision(12) << E1 << "\n";
        std::cout << "  Error energía: " << std::setprecision(4) << dE;
        if (dE < 1e-6) std::cout << " ✅\n";
        else std::cout << " ❌\n";
        
        std::cout << "  Error separación: " << dr;
        if (dr < 1e-6) std::cout << " ✅\n";
        else std::cout << " ❌\n";
        
        // Si algún paso funciona bien, el test es exitoso
        if (dE < 1e-6 && dr < 1e-6) {
            std::cout << "\n✅ KS funciona correctamente con dt = " << ks_dt << "\n";
            return true;
        }
    }
    
    std::cout << "\n❌ KS NO conserva la órbita para ningún paso probado\n";
    return false;
}

bool test_ks_phase() {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "🔬 TEST: Fase después de un período\n";
    std::cout << std::string(60, '=') << "\n";
    
    // Misma configuración
    double m1 = 1.0, m2 = 1.0;
    double a = 1.0;
    double M = m1 + m2;
    double v = std::sqrt(M / a) * 0.5;
    
    Body body1{{-a/2, 0, 0}, {0, -v, 0}, m1};
    Body body2{{ a/2, 0, 0}, {0,  v, 0}, m2};
    
    BinaryState binary(body1, body2);
    double T = kepler_period(a, M);
    
    // Posición inicial
    Vec3 r0 = binary.relative_position();
    std::cout << "Posición relativa inicial: (" 
              << r0.x << ", " << r0.y << ", " << r0.z << ")\n";
    
    KSIntegrator ks(1e-4);
    ks.integrate(binary, T);
    
    // Después de un período, debería volver a la misma posición
    Vec3 r1 = binary.relative_position();
    std::cout << "Posición relativa después de 1 período: ("
              << r1.x << ", " << r1.y << ", " << r1.z << ")\n";
    
    double error_pos = norm(r1 - r0);
    std::cout << "Error posición: " << error_pos << "\n";
    
    return error_pos < 1e-6;
}

int main() {
    std::cout << std::setprecision(12);
    
    bool test1 = test_ks_conservation();
    bool test2 = test_ks_phase();
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "RESULTADOS FINALES:\n";
    std::cout << "  Conservación: " << (test1 ? "✅" : "❌") << "\n";
    std::cout << "  Fase correcta: " << (test2 ? "✅" : "❌") << "\n";
    std::cout << std::string(60, '=') << "\n";
    
    if (test1 && test2) {
        std::cout << "🎉 ¡KS funciona correctamente!\n";
        return 0;
    } else {
        std::cout << "🔧 KS necesita correcciones\n";
        return 1;
    }
}