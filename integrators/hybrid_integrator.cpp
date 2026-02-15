// integrators/hybrid_integrator.cpp
#include "hybrid_integrator.h"

#include <algorithm>
#include <iostream>

#include "nbody_system.h"
#include "vec3.h"

// ============================================================================
// CONSTRUCTOR
// ============================================================================
HybridIntegrator::HybridIntegrator(
    std::unique_ptr<Integrator> far_integrator,
    double r_close_, //Sigue usando r_close=0.5 (KS desactivado)
    double ks_dt,
    RegimeLogger* logger_
)
: far(std::move(far_integrator))
, ks(ks_dt)
, r_close(r_close_)
, logger(logger_) {}

// ============================================================================
// DETECCIÓN FÍSICA DE UN PAR
// ============================================================================
bool HybridIntegrator::is_bound_binary(
    const NBodySystem& system,
    int i,
    int j,
    double& binding_energy
) const {
    const auto& a = system.bodies[i];
    const auto& b = system.bodies[j];

    Vec3 r = b.position - a.position;
    Vec3 v = b.velocity - a.velocity;

    double rnorm = norm(r);
    if (rnorm > r_close) return false;

    double mu = (a.mass * b.mass) / (a.mass + b.mass);
    double kinetic = 0.5 * mu * dot(v, v);
    double potential = - (a.mass * b.mass) / rnorm;

    binding_energy = kinetic + potential;
    return binding_energy < 0.0;
}

// ============================================================================
// DETECCIÓN DE TODOS LOS PARES LIGADOS
// ============================================================================
std::vector<BinaryPair>
HybridIntegrator::detect_all_bound_pairs(
    const NBodySystem& system
) const {
    std::vector<BinaryPair> all_pairs;
    const int N = static_cast<int>(system.bodies.size());

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double E;
            if (is_bound_binary(system, i, j, E)) {
                all_pairs.push_back({i, j, E});
            }
        }
    }
    return all_pairs;
}

// ============================================================================
// SELECCIÓN ÓPTIMA SIN SOLAPAMIENTOS
// ============================================================================
std::vector<BinaryPair>
HybridIntegrator::select_optimal_binaries(
    std::vector<BinaryPair>& candidates,
    std::vector<bool>& used
) const {
    // 1. Ordenar por energía más negativa (más ligados primero)
    std::sort(
        candidates.begin(),
        candidates.end(),
        [](const BinaryPair& a, const BinaryPair& b) {
            return a.binding_energy < b.binding_energy;
        }
    );

    // 2. Selección greedy óptima
    std::vector<BinaryPair> selected;
    for (const auto& pair : candidates) {
        if (!used[pair.i] && !used[pair.j]) {
            used[pair.i] = true;
            used[pair.j] = true;
            selected.push_back(pair);
        }
    }

    return selected;
}

// ============================================================================
// SEPARACIÓN DE FUERZAS (para handshake)
// ============================================================================
void HybridIntegrator::compute_external_forces(
    const NBodySystem& system,
    const std::vector<bool>& in_binary,
    std::vector<Vec3>& acc_binary,
    std::vector<Vec3>& acc_field
) const {
    const int N = static_cast<int>(system.bodies.size());
    
    acc_binary.assign(N, {0, 0, 0});
    acc_field.assign(N, {0, 0, 0});
    
    // Campo → Binarias
    for (int i = 0; i < N; ++i) {
        if (!in_binary[i]) continue;
        
        for (int j = 0; j < N; ++j) {
            if (in_binary[j]) continue;
            
            Vec3 r = system.bodies[j].position - system.bodies[i].position;
            double d = norm(r) + 1e-9;
            double factor = system.G * system.bodies[j].mass / (d * d * d);
            
            acc_binary[i] = acc_binary[i] + factor * r;
        }
    }
    
    // Binarias → Campo
    for (int i = 0; i < N; ++i) {
        if (in_binary[i]) continue;
        
        for (int j = 0; j < N; ++j) {
            if (!in_binary[j]) continue;
            
            Vec3 r = system.bodies[j].position - system.bodies[i].position;
            double d = norm(r) + 1e-9;
            double factor = system.G * system.bodies[j].mass / (d * d * d);
            
            acc_field[i] = acc_field[i] + factor * r;
        }
    }
}

// ============================================================================
// PASO PRINCIPAL DEL INTEGRADOR HÍBRIDO
// ============================================================================
void HybridIntegrator::step(
    NBodySystem& system,
    double dt,
    const std::vector<bool>& /*unused*/
) {
    const int N = static_cast<int>(system.bodies.size());
    
    // ========================================================================
    // FASE 0: PREPARACIÓN
    // ========================================================================
    std::vector<bool> in_binary(N, false);
    
    // Guardar estado inicial
    std::vector<Vec3> x0(N), v0(N);
    std::vector<double> m0(N);
    for (int i = 0; i < N; ++i) {
        x0[i] = system.bodies[i].position;
        v0[i] = system.bodies[i].velocity;
        m0[i] = system.bodies[i].mass;
    }
    
    // ========================================================================
    // FASE 1: PRIMER KICK (medio paso para TODOS)
    // ========================================================================
    auto acc_full = system.compute_accelerations();
    for (int i = 0; i < N; ++i) {
        system.bodies[i].velocity = v0[i] + 0.5 * dt * acc_full[i];
    }
    
    // ========================================================================
    // FASE 2: DETECCIÓN Y SELECCIÓN ÓPTIMA DE BINARIAS
    // ========================================================================
    // 2.1 Detectar TODOS los pares ligados
    auto all_pairs = detect_all_bound_pairs(system);
    
    // 2.2 Seleccionar los mejores sin solapamientos
    auto selected_binaries = select_optimal_binaries(all_pairs, in_binary);
    
    // ========================================================================
    // FASE 3: INTEGRACIÓN DE BINARIAS CON KS
    // ========================================================================
    for (const auto& bin : selected_binaries) {
        if (logger) {
            logger->log({
                static_cast<int>(step_counter),
                bin.i,
                bin.j,
                "ENTER_KS"
            });
        }
        
        // Crear y integrar binaria
        BinaryState state(system.bodies[bin.i], system.bodies[bin.j]);
        ks.integrate(state, dt);  // KS subcicla internamente
        state.write_back(system.bodies[bin.i], system.bodies[bin.j]);
        
        if (logger) {
            logger->log({
                static_cast<int>(step_counter),
                bin.i,
                bin.j,
                "EXIT_KS"
            });
        }
    }
    
    // ========================================================================
    // FASE 4: DRIFT DE CUERPOS DEL CAMPO
    // ========================================================================
    for (int i = 0; i < N; ++i) {
        if (!in_binary[i]) {
            system.bodies[i].position = x0[i] + dt * system.bodies[i].velocity;
        }
    }
    
    // ========================================================================
    // FASE 5: SEGUNDO KICK
    // ========================================================================
    system.invalidate_accelerations();
    auto acc_new = system.compute_accelerations();
    
    for (int i = 0; i < N; ++i) {
        if (!in_binary[i]) {
            system.bodies[i].velocity = v0[i] + 0.5 * dt * (acc_full[i] + acc_new[i]);
        }
    }
    
    // ========================================================================
    // FASE 6: VERIFICACIÓN OPCIONAL
    // ========================================================================
    #ifdef DEBUG_HYBRID
    double E_before = 0;
    for (int i = 0; i < N; ++i) {
        E_before += 0.5 * m0[i] * dot(v0[i], v0[i]);
        for (int j = i+1; j < N; ++j) {
            Vec3 r0 = x0[j] - x0[i];
            E_before -= system.G * m0[i] * m0[j] / (norm(r0) + 1e-9);
        }
    }
    
    double E_after = system.total_energy();
    double dE = std::abs(E_after - E_before) / std::abs(E_before);
    if (dE > 1e-6) {
        std::cout << "Warning: Paso " << step_counter 
                  << ", dE/E = " << dE << std::endl;
    }
    #endif
    
    ++step_counter;
}