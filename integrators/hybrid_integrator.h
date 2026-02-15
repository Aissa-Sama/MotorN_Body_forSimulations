// integrators/hybrid_integrator.h
#pragma once

#include <memory>
#include <vector>
#include "integrator.h"
#include "binary_state.h"
#include "ks_integrator.h"
#include "regime_logger.h"

// =========================
// Par binario candidato
// =========================
struct BinaryPair {
    int i;
    int j;
    double binding_energy;  // más negativo = más ligado
};

class HybridIntegrator : public Integrator {
public:
    HybridIntegrator(
        std::unique_ptr<Integrator> far_integrator,
        double r_close,
        double ks_internal_dt,
        RegimeLogger* logger = nullptr
    );

    void step(
        NBodySystem& system,
        double dt,
        const std::vector<bool>& used
    ) override;

private:
    // Detección física de un par
    bool is_bound_binary(
        const NBodySystem& system,
        int i,
        int j,
        double& binding_energy
    ) const;
    
    // Detecta TODOS los pares ligados (sin filtro)
    std::vector<BinaryPair> detect_all_bound_pairs(
        const NBodySystem& system
    ) const;
    
    // Selecciona conjunto óptimo sin solapamientos
    std::vector<BinaryPair> select_optimal_binaries(
        std::vector<BinaryPair>& candidates,
        std::vector<bool>& used
    ) const;
    
    // Separación de fuerzas (para handshake)
    void compute_external_forces(
        const NBodySystem& system,
        const std::vector<bool>& in_binary,
        std::vector<Vec3>& acc_binary,
        std::vector<Vec3>& acc_field
    ) const;

    // Miembros
    std::unique_ptr<Integrator> far;
    KSIntegrator ks;
    double r_close;
    RegimeLogger* logger;
    std::size_t step_counter = 0;
};

// HybridIntegrator combina dos métodos de integración: Velocity Verlet para la evolución global y RK45 para encuentros cercanos.
// El método step() decide qué integrador usar en cada paso, dependiendo de si hay un encuentro cercano entre cuerpos.

//  main
//    ↓
//  HybridIntegrator
//    ↓
//  ┌───────────────┬──────────────────┐
//  │ régimen global│ régimen local    │
//  │ verlet        │ RK45 / KS        │
//  │               │                  │
//  │ conserva H    │ resuelve         │
//  │               │ encuentros       │
//  └───────────────┴──────────────────┘
//    ↓
//  NBodySystem actualizado
