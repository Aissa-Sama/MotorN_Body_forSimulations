// tests/test_chain3_state.cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <vector>
#include "nbody_system.h"
#include "chain3_integrator.h"
#include "initial_conditions.h"

/**
 * @brief Test de las transformaciones físico ↔ cadena.
 * 
 * Verifica que un sistema de 3 cuerpos pueda convertirse a variables de cadena
 * y luego reconstruirse sin pérdida de precisión.
 * La comparación se basa en que exista una permutación de los cuerpos
 * que haga coincidir posiciones y velocidades con el sistema original.
 */
int main() {
    std::cout << std::setprecision(15);
    std::cout << "=== Test de transformaciones Físico ↔ Cadena ===\n";

    // Crear un sistema de 3 cuerpos (figura-8)
    NBodySystem system = InitialConditions::figure_eight();
    
    std::cout << "\nSistema original:\n";
    for (size_t i = 0; i < system.bodies.size(); ++i) {
        const auto& b = system.bodies[i];
        std::cout << "Cuerpo " << i << ": masa=" << b.mass
                  << "  pos=(" << b.position.x << ", " << b.position.y << ", " << b.position.z << ")"
                  << "  vel=(" << b.velocity.x << ", " << b.velocity.y << ", " << b.velocity.z << ")\n";
    }

    // Inicializar el integrador de cadena
    Chain3Integrator chain3(1e-4);

    // Convertir a variables de cadena (usando el orden dado: 0,1,2)
    Chain3State state = chain3.initialize(system, 0, 1, 2);

    std::cout << "\nVariables de cadena:\n";
    std::cout << "Q1 = (" << state.Q1.x << ", " << state.Q1.y << ", " << state.Q1.z << ", " << state.Q1.w << ")\n";
    std::cout << "P1 = (" << state.P1.x << ", " << state.P1.y << ", " << state.P1.z << ", " << state.P1.w << ")\n";
    std::cout << "Q2 = (" << state.Q2.x << ", " << state.Q2.y << ", " << state.Q2.z << ", " << state.Q2.w << ")\n";
    std::cout << "P2 = (" << state.P2.x << ", " << state.P2.y << ", " << state.P2.z << ", " << state.P2.w << ")\n";
    std::cout << "Energía interna = " << state.energy << "\n";
    std::cout << "CM: pos=(" << state.cm_pos.x << ", " << state.cm_pos.y << ", " << state.cm_pos.z << ")\n";
    std::cout << "CM: vel=(" << state.cm_vel.x << ", " << state.cm_vel.y << ", " << state.cm_vel.z << ")\n";

    // Reconstruir el sistema físico
    NBodySystem system2;
    system2.bodies.resize(3);
    system2.bodies[0].mass = state.m1();
    system2.bodies[1].mass = state.m2();
    system2.bodies[2].mass = state.m3();
    chain3.write_back(state, system2, 0, 1, 2);

    std::cout << "\nSistema reconstruido:\n";
    for (size_t i = 0; i < system2.bodies.size(); ++i) {
        const auto& b = system2.bodies[i];
        std::cout << "Cuerpo " << i << ": masa=" << b.mass
                  << "  pos=(" << b.position.x << ", " << b.position.y << ", " << b.position.z << ")"
                  << "  vel=(" << b.velocity.x << ", " << b.velocity.y << ", " << b.velocity.z << ")\n";
    }

    // ========================================================================
    // VERIFICACIÓN POR PERMUTACIONES
    // ========================================================================

    // Guardar posiciones originales
    std::vector<Vec3> orig_pos = {
        system.bodies[0].position,
        system.bodies[1].position,
        system.bodies[2].position
    };

    // Guardar posiciones reconstruidas
    std::vector<Vec3> rec_pos = {
        system2.bodies[0].position,
        system2.bodies[1].position,
        system2.bodies[2].position
    };

    // Probar las 6 permutaciones
    int best_perm[3] = {0,1,2};
    double best_error = 1e99;
    int perms[6][3] = {
        {0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}
    };

    for (int p = 0; p < 6; ++p) {
        double error = 0.0;
        for (int i = 0; i < 3; ++i) {
            error += (rec_pos[perms[p][i]] - orig_pos[i]).norm();
        }
        if (error < best_error) {
            best_error = error;
            for (int i = 0; i < 3; ++i) best_perm[i] = perms[p][i];
        }
    }

    std::cout << "\n==================================================\n";
    std::cout << "VERIFICACIÓN POR PERMUTACIONES:\n";
    std::cout << "==================================================\n";
    std::cout << "Mejor permutación encontrada: ["
              << best_perm[0] << ", " << best_perm[1] << ", " << best_perm[2] << "]\n";
    std::cout << "Error de posición tras permutación: " << best_error << "\n";

    // Aplicar la mejor permutación a las posiciones y velocidades para las siguientes comprobaciones
    NBodySystem system2_permuted;
    system2_permuted.bodies.resize(3);
    for (int i = 0; i < 3; ++i) {
        system2_permuted.bodies[i] = system2.bodies[best_perm[i]];
    }

    // ========================================================================
    // VERIFICACIÓN DE INVARIANTES FÍSICAS (con la permutación aplicada)
    // ========================================================================

    // 1. Centro de masa y momento total (deben conservarse exactamente)
    double M = system.bodies[0].mass + system.bodies[1].mass + system.bodies[2].mass;
    Vec3 cm_pos_orig = (system.bodies[0].mass * system.bodies[0].position +
                        system.bodies[1].mass * system.bodies[1].position +
                        system.bodies[2].mass * system.bodies[2].position) / M;
    Vec3 cm_vel_orig = (system.bodies[0].mass * system.bodies[0].velocity +
                        system.bodies[1].mass * system.bodies[1].velocity +
                        system.bodies[2].mass * system.bodies[2].velocity) / M;

    Vec3 cm_pos_rec = (system2_permuted.bodies[0].mass * system2_permuted.bodies[0].position +
                       system2_permuted.bodies[1].mass * system2_permuted.bodies[1].position +
                       system2_permuted.bodies[2].mass * system2_permuted.bodies[2].position) / M;
    Vec3 cm_vel_rec = (system2_permuted.bodies[0].mass * system2_permuted.bodies[0].velocity +
                       system2_permuted.bodies[1].mass * system2_permuted.bodies[1].velocity +
                       system2_permuted.bodies[2].mass * system2_permuted.bodies[2].velocity) / M;

    double cm_pos_error = (cm_pos_rec - cm_pos_orig).norm();
    double cm_vel_error = (cm_vel_rec - cm_vel_orig).norm();

    // 2. Energía total (debe conservarse)
    double E_orig = system.total_energy();
    double E_rec = system2_permuted.total_energy();
    double energy_error = std::abs(E_rec - E_orig) / std::abs(E_orig);

    // 3. Todas las distancias entre pares (sin importar el orden)
    std::vector<double> dist_orig;
    for (int i = 0; i < 3; ++i) {
        for (int j = i+1; j < 3; ++j) {
            dist_orig.push_back((system.bodies[j].position - system.bodies[i].position).norm());
        }
    }
    std::sort(dist_orig.begin(), dist_orig.end());

    std::vector<double> dist_rec;
    for (int i = 0; i < 3; ++i) {
        for (int j = i+1; j < 3; ++j) {
            dist_rec.push_back((system2_permuted.bodies[j].position - system2_permuted.bodies[i].position).norm());
        }
    }
    std::sort(dist_rec.begin(), dist_rec.end());

    double max_dist_error = 0.0;
    for (size_t i = 0; i < dist_orig.size(); ++i) {
        max_dist_error = std::max(max_dist_error, std::abs(dist_rec[i] - dist_orig[i]));
    }

    // 4. Todas las velocidades relativas (normas)
    std::vector<double> vrel_orig;
    for (int i = 0; i < 3; ++i) {
        for (int j = i+1; j < 3; ++j) {
            vrel_orig.push_back((system.bodies[j].velocity - system.bodies[i].velocity).norm());
        }
    }
    std::sort(vrel_orig.begin(), vrel_orig.end());

    std::vector<double> vrel_rec;
    for (int i = 0; i < 3; ++i) {
        for (int j = i+1; j < 3; ++j) {
            vrel_rec.push_back((system2_permuted.bodies[j].velocity - system2_permuted.bodies[i].velocity).norm());
        }
    }
    std::sort(vrel_rec.begin(), vrel_rec.end());

    double max_vrel_error = 0.0;
    for (size_t i = 0; i < vrel_orig.size(); ++i) {
        max_vrel_error = std::max(max_vrel_error, std::abs(vrel_rec[i] - vrel_orig[i]));
    }

    // 5. Momento angular total (vector)
    Vec3 L_orig = system.total_angular_momentum();
    Vec3 L_rec = system2_permuted.total_angular_momentum();
    double L_error = (L_rec - L_orig).norm();

    // ========================================================================
    // REPORTE DE RESULTADOS
    // ========================================================================

    std::cout << "\n" << std::string(50, '=') << "\n";
    std::cout << "RESULTADOS DE LA VERIFICACIÓN (con permutación):\n";
    std::cout << std::string(50, '=') << "\n";

    std::cout << "Error en posición del CM: " << cm_pos_error << "\n";
    std::cout << "Error en velocidad del CM: " << cm_vel_error << "\n";
    std::cout << "Error máximo en distancias (ordenadas): " << max_dist_error << "\n";
    std::cout << "Error máximo en velocidades relativas (ordenadas): " << max_vrel_error << "\n";
    std::cout << "Error en momento angular total: " << L_error << "\n";
    std::cout << "Error relativo en energía total: " << energy_error << "\n";

    std::cout << std::string(50, '-') << "\n";

    // Criterios de éxito (tolerancia generosa para la Fase 2)
    bool cm_ok = (cm_pos_error < 1e-12 && cm_vel_error < 1e-12);
    bool dist_ok = (max_dist_error < 1e-12);
    bool vrel_ok = (max_vrel_error < 1e-12);
    bool L_ok = (L_error < 1e-12);
    bool energy_ok = (energy_error < 1e-12);

    if (cm_ok && dist_ok && vrel_ok && L_ok && energy_ok) {
        std::cout << "✅ Test pasado: todas las invariantes se conservan.\n";
        return 0;
    } else {
        std::cout << "❌ Test fallado: algunas invariantes no se conservan.\n";
        if (!cm_ok) std::cout << "   - Centro de masa no se conserva.\n";
        if (!dist_ok) std::cout << "   - Distancias no se conservan.\n";
        if (!vrel_ok) std::cout << "   - Velocidades relativas no se conservan.\n";
        if (!L_ok) std::cout << "   - Momento angular no se conserva.\n";
        if (!energy_ok) std::cout << "   - Energía total no se conserva.\n";
        return 1;
    }
}