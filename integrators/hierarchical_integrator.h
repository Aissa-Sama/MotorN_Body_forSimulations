// integrators/hierarchical_integrator.h
#pragma once
#include <memory>
#include <vector>
#include <map>
#include <tuple>
#include <functional>
#include "integrator.h"
#include "hierarchy_node.h"
#include "hierarchy_builder.h"
#include "chain3_integrator.h"
#include "archain3_integrator.h"
#include "archain3_state.h"
#include "archain_n_ks_integrator.h"
#include "ks_perturbed_integrator.h"
#include "ks_integrator.h"
#include "regime_logger.h"
#include "nbody_system.h"
#include "block_timestep.h"    // Fase 7B
#include "merger_engine.h"     // Fase 7C

// ============================================================================
// HIERARCHICAL INTEGRATOR — Ruta B
//
// Cada paso:
//   0. [Fase 7C] MergerEngine detecta y ejecuta fusiones/TDEs (opcional)
//   1. HierarchyBuilder::build() → árbol del instante actual
//   2. Recorrer el árbol e integrar cada nodo con su método especializado:
//        LEAF             → far (Leapfrog) [con block timestep si activo]
//        PAIR_KS          → KS simple o KS perturbado (tidal_parameter)
//        TRIPLE_CHAIN     → Chain3Integrator  (sep ≥ ar_chain_threshold)
//        GROUP_AR_CHAIN   → ARChainNKSIntegrator (sep < ar_chain_threshold)
//        COMPOSITE        → integrar hijos recursivamente
//
// ── FASES ACUMULADAS ─────────────────────────────────────────────────────────
//
// Fase 6A: HierarchyBuilder marca pn_active; despacho a ARChainNPNBSIntegrator.
// Fase 6B: ARChainNKSIntegrator con paso ds_ks = η·√(sep·r_ks)/Ω.
// Fase 7A: GROUP_AR_CHAIN para N≥3 (fusión de TRIPLE_AR_CHAIN).
// Fase 7B: block_timestep_ — dt individual por cuerpo LEAF (Aarseth 2003).
// Fase 7C: merger_engine_ — fusiones físicas, TDE, captura Schwarzschild.
//
// ── MODOS DE OPERACIÓN ───────────────────────────────────────────────────────
//
// 1. step(system, dt) — modo bloque: avanza dt. Uso general, N > 3 o mixto.
// 2. step_to(system, t_final) — modo directo: sin cortes externos.
//    Delega a ARChainNKSIntegrator::integrate_to() si el árbol es un único
//    GROUP_AR_CHAIN. Más preciso para sistemas totalmente acoplados.
//
// ── PERSISTENCIA DEL ESTADO AR-CHAIN ─────────────────────────────────────────
// ar_chain_n_states_: map<string, ARChainNState> indexado por clave canónica
// (índices ordenados, concatenados como "i_j_k_...").
// El estado persiste entre pasos para sobrevivir reconstrucciones del árbol.
// ============================================================================

class HierarchicalIntegrator : public Integrator {
public:

    HierarchicalIntegrator(
        std::unique_ptr<Integrator>     far_integrator,
        double                          r_ks_threshold,
        double                          ks_internal_dt,
        const HierarchyBuilder::Params& builder_params  = HierarchyBuilder::Params{},
        RegimeLogger*                   logger           = nullptr,
        // Fase 7B — block timestep
        bool                            enable_block_ts  = false,
        const BlockTimestep::Params&    bt_params        = BlockTimestep::Params{},
        // Fase 7C — merger engine
        bool                            enable_mergers   = false,
        const tidal::MergerEngine::Params& merger_params = tidal::MergerEngine::Params{}
    );

    // ── Interfaz principal ───────────────────────────────────────────────────

    /** Modo bloque: avanza dt. Uso general. */
    void step(
        NBodySystem&              system,
        double                    dt,
        const std::vector<bool>&  used
    ) override;

    /**
     * Modo directo: integra de t=0 a t_final sin cortes externos.
     * Si el árbol entero es un único GROUP_AR_CHAIN, delega directamente
     * a ARChainNKSIntegrator::integrate_to().
     */
    void step_to(
        NBodySystem&                    system,
        double                          t_final,
        std::function<void(double)>     logger_cb = nullptr
    );

    // ── Acceso al árbol del último paso ─────────────────────────────────────
    const HierarchyNode* last_tree() const { return last_root.get(); }

    // ── Acceso al merger engine (para configurar callbacks externos) ─────────
    tidal::MergerEngine& merger_engine() { return merger_engine_; }

private:

    // ── Detección del modo directo ───────────────────────────────────────────
    // Devuelve true si el árbol raíz es un único GROUP_AR_CHAIN con índices
    // (out_indices). También acepta nodos TRIPLE_AR_CHAIN (alias).
    bool is_pure_group_ar_chain(const HierarchyNode&  root,
                                std::vector<int>&      out_indices) const;

    // ── Despachador recursivo ────────────────────────────────────────────────
    void integrate_node(
        HierarchyNode&        node,
        NBodySystem&          system,
        double                dt,
        std::vector<bool>&    in_subsystem
    );

    void integrate_leaf(
        HierarchyNode&        node,
        NBodySystem&          system,
        double                dt,
        std::vector<bool>&    in_subsystem
    );

    void integrate_pair_ks(
        HierarchyNode&  node,
        NBodySystem&    system,
        double          dt
    );

    void integrate_triple_chain(
        HierarchyNode&  node,
        NBodySystem&    system,
        double          dt
    );

    void integrate_group_ar_chain(
        HierarchyNode&  node,
        NBodySystem&    system,
        double          dt
    );

    void integrate_composite(
        HierarchyNode&      node,
        NBodySystem&        system,
        double              dt,
        std::vector<bool>&  in_subsystem
    );

    // ── Clave canónica para caché AR-chain ───────────────────────────────────
    // Índices ordenados concatenados: "2_5_7" para {7,2,5}.
    static std::string make_key(std::vector<int> indices) {
        std::sort(indices.begin(), indices.end());
        std::string key;
        for (int i : indices) {
            if (!key.empty()) key += '_';
            key += std::to_string(i);
        }
        return key;
    }

    // ── Miembros ─────────────────────────────────────────────────────────────
    std::unique_ptr<Integrator>   far;
    KSIntegrator                  ks_simple;
    KSPerturbedIntegrator         ks_perturbed;
    Chain3Integrator              chain3;
    ARChainNKSIntegrator          ar_chain_ks_;   // Fases 6B + 7A
    HierarchyBuilder              builder;
    double                        tidal_threshold;
    RegimeLogger*                 logger;
    std::size_t                   step_counter = 0;
    double                        current_time_ = 0.0;  // tiempo físico acumulado

    std::unique_ptr<HierarchyNode>           last_root;
    std::map<std::string, ARChainNState>     ar_chain_n_states_;

    // Fase 7B — block timestep
    bool               enable_block_ts_;
    BlockTimestep      block_ts_;

    // Fase 7C — merger engine
    bool               enable_mergers_;
    tidal::MergerEngine merger_engine_;
};
