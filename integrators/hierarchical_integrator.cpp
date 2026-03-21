// integrators/hierarchical_integrator.cpp
#include "hierarchical_integrator.h"
#include "binary_state.h"
#include "leapfrog_integrator.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <functional>
#include <numeric>

// ============================================================================
// CONSTRUCTOR
// ============================================================================
HierarchicalIntegrator::HierarchicalIntegrator(
    std::unique_ptr<Integrator>         far_integrator,
    double                              r_ks_threshold,
    double                              ks_internal_dt,
    const HierarchyBuilder::Params&     builder_params,
    RegimeLogger*                       logger_,
    bool                                enable_block_ts,
    const BlockTimestep::Params&        bt_params,
    bool                                enable_mergers,
    const tidal::MergerEngine::Params&  merger_params
)
    : far(std::move(far_integrator))
    , ks_simple(ks_internal_dt)
    , ks_perturbed(ks_internal_dt)
    , chain3(1e-4)
    , ar_chain_ks_(builder_params.ar_chain_eta)
    , builder(builder_params)
    , tidal_threshold(builder_params.tidal_threshold)
    , logger(logger_)
    , enable_block_ts_(enable_block_ts)
    , block_ts_(bt_params)
    , enable_mergers_(enable_mergers)
    , merger_engine_(merger_params)
{
    (void)r_ks_threshold;
}

// ============================================================================
// IS_PURE_GROUP_AR_CHAIN
// Detecta si el árbol raíz es un único nodo GROUP_AR_CHAIN (o TRIPLE_AR_CHAIN,
// que es alias). Rellena out_indices con los índices del grupo.
// ============================================================================
bool HierarchicalIntegrator::is_pure_group_ar_chain(
    const HierarchyNode& root,
    std::vector<int>&    out_indices) const
{
    if (root.type != HierarchyNode::Type::GROUP_AR_CHAIN &&
        root.type != HierarchyNode::Type::TRIPLE_AR_CHAIN)
        return false;
    if (root.body_indices.empty()) return false;

    out_indices = root.body_indices;
    return true;
}

// ============================================================================
// STEP — modo bloque (uso general)
// ============================================================================
void HierarchicalIntegrator::step(
    NBodySystem&              system,
    double                    dt,
    const std::vector<bool>&  /*external_mask*/
)
{
    const int N = static_cast<int>(system.bodies.size());

    // ── Fase 7C: detectar y ejecutar fusiones ANTES de construir el árbol ───
    if (enable_mergers_) {
        int n_ev = merger_engine_.process(system, current_time_);
        if (n_ev > 0) system.invalidate_accelerations();
    }

    // ── Construir árbol de jerarquía ─────────────────────────────────────────
    last_root = builder.build(system);

    std::vector<bool> in_subsystem(system.bodies.size(), false);
    integrate_node(*last_root, system, dt, in_subsystem);

    // ── Fase 7B: block timestep para LEAFs (o leapfrog uniforme) ────────────
    if (enable_block_ts_) {
        // Calcular aceleraciones y jerks para los cuerpos libres
        auto accs  = system.compute_accelerations();
        auto jerks = block_ts_.compute_jerks(system, accs);
        auto dts   = block_ts_.assign(system, accs, jerks, in_subsystem, dt);

        // step_block: DKD con dt individual
        if (auto* lf = dynamic_cast<LeapfrogIntegrator*>(far.get())) {
            lf->step_block(system, dts, in_subsystem);
        } else {
            // Fallback: leapfrog uniforme con dt
            far->step(system, dt, in_subsystem);
        }
    } else {
        far->step(system, dt, in_subsystem);
    }

    current_time_ += dt;
    ++step_counter;
}

// ============================================================================
// STEP_TO — modo directo (sistema completo, sin cortes externos)
// ============================================================================
void HierarchicalIntegrator::step_to(
    NBodySystem&                    system,
    double                          t_final,
    std::function<void(double)>     logger_cb
)
{
    // ── Fase 7C: fusiones antes de comenzar la integración ──────────────────
    if (enable_mergers_) {
        int n_ev = merger_engine_.process(system, current_time_);
        if (n_ev > 0) system.invalidate_accelerations();
    }

    // ── Construir árbol y detectar si es un único grupo AR-chain ─────────────
    last_root = builder.build(system);

    std::vector<int> group_indices;
    if (!is_pure_group_ar_chain(*last_root, group_indices)) {
        // Sistema mixto — caer al modo step() con dt heurístico
        const double dt_fallback = t_final / 1000.0;
        double t_cur = current_time_;
        while (t_cur < current_time_ + t_final - 1e-14) {
            const double dt = std::min(dt_fallback,
                                       current_time_ + t_final - t_cur);
            step(system, dt, {});
            t_cur += dt;
            if (logger_cb) logger_cb(t_cur);
        }
        return;
    }

    // ── Integración directa: GROUP_AR_CHAIN sin cortes externos ─────────────
    if (logger)
        logger->log({static_cast<int>(step_counter),
                     group_indices[0], group_indices.back(),
                     "ENTER_AR_CHAIN_DIRECT"});

    ARChainNState state = ar_chain_ks_.initialize(system, group_indices);

    // Clave canónica para el caché
    std::string key = make_key(group_indices);

    // Integrar directamente hasta t_final (tiempo absoluto desde t=0 del estado)
    ar_chain_ks_.integrate_to(state, t_final);

    // Escribir resultado al sistema físico
    ar_chain_ks_.write_back(state, system, group_indices);
    ar_chain_n_states_[key] = state;

    if (logger)
        logger->log({static_cast<int>(step_counter),
                     group_indices[0], group_indices.back(),
                     "EXIT_AR_CHAIN_DIRECT"});

    current_time_ += t_final;
    ++step_counter;

    if (logger_cb) logger_cb(current_time_);
}

// ============================================================================
// DESPACHADOR RECURSIVO
// ============================================================================
void HierarchicalIntegrator::integrate_node(
    HierarchyNode&      node,
    NBodySystem&        system,
    double              dt,
    std::vector<bool>&  in_subsystem
)
{
    switch (node.type) {
        case HierarchyNode::Type::LEAF:
            integrate_leaf(node, system, dt, in_subsystem);
            break;
        case HierarchyNode::Type::PAIR_KS:
            integrate_pair_ks(node, system, dt);
            for (int idx : node.body_indices) in_subsystem[idx] = true;
            break;
        case HierarchyNode::Type::TRIPLE_CHAIN:
            integrate_triple_chain(node, system, dt);
            for (int idx : node.body_indices) in_subsystem[idx] = true;
            break;
        case HierarchyNode::Type::GROUP_AR_CHAIN:
        // TRIPLE_AR_CHAIN es alias de GROUP_AR_CHAIN en el enum — mismo case
            integrate_group_ar_chain(node, system, dt);
            for (int idx : node.body_indices) in_subsystem[idx] = true;
            break;
        case HierarchyNode::Type::COMPOSITE:
            integrate_composite(node, system, dt, in_subsystem);
            break;
        default:
            break;
    }
}

// ============================================================================
// LEAF: marcar como libre (far lo integrará)
// ============================================================================
void HierarchicalIntegrator::integrate_leaf(
    HierarchyNode&      /*node*/,
    NBodySystem&        /*system*/,
    double              /*dt*/,
    std::vector<bool>&  /*in_subsystem*/
)
{
    // Los LEAFs son procesados por far->step() al final de step().
    // No se marcan en in_subsystem aquí.
}

// ============================================================================
// PAIR_KS: KS simple o perturbado según tidal_parameter
// ============================================================================
void HierarchicalIntegrator::integrate_pair_ks(
    HierarchyNode&  node,
    NBodySystem&    system,
    double          dt
)
{
    const int i = node.body_indices[0];
    const int j = node.body_indices[1];

    if (logger) {
        logger->log({
            static_cast<int>(step_counter), i, j,
            node.tidal_parameter < tidal_threshold
                ? "ENTER_KS_SIMPLE" : "ENTER_KS_PERTURBED"
        });
    }

    BinaryState state(system.bodies[i], system.bodies[j]);

    if (node.tidal_parameter < tidal_threshold)
        ks_simple.integrate(state, dt);
    else
        ks_perturbed.integrate_perturbed(state, dt, system, i, j);

    state.write_back(system.bodies[i], system.bodies[j]);
}

// ============================================================================
// TRIPLE_CHAIN: Chain3Integrator (KS-chain, separación moderada)
// ============================================================================
void HierarchicalIntegrator::integrate_triple_chain(
    HierarchyNode&  node,
    NBodySystem&    system,
    double          dt
)
{
    const int i = node.body_indices[0];
    const int j = node.body_indices[1];
    const int k = node.body_indices[2];

    if (logger)
        logger->log({static_cast<int>(step_counter), i, j, "ENTER_CHAIN3"});

    Chain3State state = chain3.initialize(system, i, j, k);

    double t_target   = state.cm_time + dt;
    double t_achieved = 0.0;
    IntegrationParams params;
    params.abs_tol  = 1e-10;
    params.min_dtau = 1e-8;
    params.max_dtau = 1e-1;

    chain3.integrate(state, t_target, t_achieved, params, system);
    chain3.write_back(state, system, i, j, k);

    if (logger)
        logger->log({static_cast<int>(step_counter), i, j, "EXIT_CHAIN3"});
}

// ============================================================================
// GROUP_AR_CHAIN: ARChainNKSIntegrator en modo bloque
//
// Modo bloque: llamado desde step() con un dt fijo.
// Para el modo directo (sin cortes), usar step_to().
//
// PERSISTENCIA DEL ESTADO:
//   El árbol se reconstruye en cada paso — los nodos son efímeros.
//   El estado ARChainNState se guarda en ar_chain_n_states_ indexado
//   por clave canónica (índices ordenados, string "i_j_k_...").
//
// MANEJO DEL CM ENTRE BLOQUES:
//   Al retomar el estado previo:
//     1. Absorber desplazamiento del CM: cm_pos += cm_vel * t_phys
//     2. Resetear el reloj interno:      t_phys  = 0
// ============================================================================
void HierarchicalIntegrator::integrate_group_ar_chain(
    HierarchyNode&  node,
    NBodySystem&    system,
    double          dt
)
{
    const std::vector<int>& indices = node.body_indices;
    if (indices.empty()) return;

    if (logger)
        logger->log({static_cast<int>(step_counter),
                     indices[0], indices.back(), "ENTER_AR_CHAIN"});

    std::string key = make_key(indices);

    ARChainNState state;
    auto it = ar_chain_n_states_.find(key);
    if (it != ar_chain_n_states_.end()) {
        state = it->second;
        // Absorber desplazamiento del CM del bloque anterior y resetear reloj
        state.cm_pos = state.cm_pos + state.cm_vel * state.t_phys;
        state.t_phys = 0.0;
    } else {
        state = ar_chain_ks_.initialize(system, indices);
    }

    ar_chain_ks_.integrate(state, dt);

    ar_chain_n_states_[key] = state;
    ar_chain_ks_.write_back(state, system, indices);

    if (logger)
        logger->log({static_cast<int>(step_counter),
                     indices[0], indices.back(), "EXIT_AR_CHAIN"});
}

// ============================================================================
// COMPOSITE: integrar hijos recursivamente
// ============================================================================
void HierarchicalIntegrator::integrate_composite(
    HierarchyNode&      node,
    NBodySystem&        system,
    double              dt,
    std::vector<bool>&  in_subsystem
)
{
    for (auto& child : node.children)
        integrate_node(*child, system, dt, in_subsystem);
}
