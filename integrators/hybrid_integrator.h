// integrators/hybrid_integrator.h
#pragma once
#include <memory>
#include <vector>
#include "integrator.h"
#include "binary_state.h"
#include "ks_perturbed_integrator.h"
#include "chain3_integrator.h"
#include "regime_logger.h"
#include "nbody_system.h"

// ============================================================================
// Par binario candidato
// ============================================================================
struct BinaryPair {
    int i;
    int j;
    double binding_energy;  // más negativo = más ligado
};

// ============================================================================
// Triple candidato (para Chain3)
// ============================================================================
struct TripleCandidate {
    int i, j, k;           // índices en NBodySystem (ordenados para la cadena)
    double binding_energy;  // energía de ligadura total del triple
};

// ============================================================================
// HybridIntegrator
//
// Arquitectura para Ruta B (HierarchyBuilder):
//   - detect_triple()         → será reemplazado por HierarchyBuilder
//   - handle_triple()         → será reemplazado por un nodo TRIPLE_CHAIN
//   - handle_binary_ks()      → será reemplazado por un nodo PAIR_KS
//   - La lógica de step() se mantiene igual; solo cambia quién detecta/clasifica
// ============================================================================
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
    // -------------------------------------------------------------------------
    // DETECCIÓN — Ruta B reemplazará estos métodos con HierarchyBuilder
    // -------------------------------------------------------------------------

    // Detecta un triple fuertemente acoplado en el sistema.
    // Devuelve true y rellena 'triple' si encuentra uno.
    // HOOK para Ruta B: HierarchyBuilder::find_triples() reemplaza esto.
    bool detect_triple(
        const NBodySystem& system,
        TripleCandidate& triple
    ) const;

    // Verifica si un par específico está ligado
    bool is_bound_binary(
        const NBodySystem& system,
        int i, int j,
        double& binding_energy
    ) const;

    // Detecta todos los pares ligados
    std::vector<BinaryPair> detect_all_bound_pairs(
        const NBodySystem& system
    ) const;

    // Selecciona conjunto óptimo sin solapamientos
    std::vector<BinaryPair> select_optimal_binaries(
        std::vector<BinaryPair>& candidates,
        std::vector<bool>& used,
        const NBodySystem& system
    ) const;

    // -------------------------------------------------------------------------
    // INTEGRACIÓN — Ruta B reemplazará estos métodos con nodos del árbol
    // -------------------------------------------------------------------------

    // Integra un triple con Chain3Integrator.
    // HOOK para Ruta B: HierarchyNode::TRIPLE_CHAIN::integrate() reemplaza esto.
    void handle_triple(
        NBodySystem& system,
        const TripleCandidate& triple,
        double dt
    );

    // Integra una binaria con KS perturbado.
    // HOOK para Ruta B: HierarchyNode::PAIR_KS::integrate() reemplaza esto.
    void handle_binary_ks(
        NBodySystem& system,
        const BinaryPair& bin,
        double dt
    );

    // Separación de fuerzas externas
    void compute_external_forces(
        const NBodySystem& system,
        const std::vector<bool>& in_binary,
        std::vector<Vec3>& acc_binary,
        std::vector<Vec3>& acc_field
    ) const;

    // -------------------------------------------------------------------------
    // CRITERIOS DE CLASIFICACIÓN
    // -------------------------------------------------------------------------

    // Verifica si tres cuerpos están mutuamente dentro de r_close
    bool is_close_triple(
        const NBodySystem& system,
        int i, int j, int k
    ) const;

    // Detecta si una binaria KS tiene un tercer cuerpo demasiado cerca
    // → debe escalarse a Chain3
    bool third_body_too_close(
        const NBodySystem& system,
        int bi, int bj,     // la binaria
        int& third_idx      // salida: índice del tercer cuerpo
    ) const;

    // -------------------------------------------------------------------------
    // MIEMBROS
    // -------------------------------------------------------------------------
    std::unique_ptr<Integrator> far;
    KSPerturbedIntegrator       ks_perturbed;
    Chain3Integrator            chain3;
    double                      r_close;
    RegimeLogger*               logger;
    std::size_t                 step_counter = 0;
};
