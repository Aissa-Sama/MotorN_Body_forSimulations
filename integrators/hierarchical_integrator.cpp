// integrators/hierarchical_integrator.cpp
#include "hierarchical_integrator.h"
#include "binary_state.h"
#include <cmath>
#include <stdexcept>

// ============================================================================
// CONSTRUCTOR
// ============================================================================
HierarchicalIntegrator::HierarchicalIntegrator(
    std::unique_ptr<Integrator> far_integrator,
    double r_ks_threshold,
    double ks_internal_dt,
    const HierarchyBuilder::Params& builder_params,
    RegimeLogger* logger_
)
    : far(std::move(far_integrator))
    , ks_simple(ks_internal_dt)
    , ks_perturbed(ks_internal_dt)
    , chain3(1e-4)
    , builder(builder_params)
    , tidal_threshold(builder_params.tidal_threshold)
    , logger(logger_)
{
    (void)r_ks_threshold;  // queda en builder_params.r_ks_threshold
}

// ============================================================================
// PASO PRINCIPAL
// ============================================================================
void HierarchicalIntegrator::step(
    NBodySystem& system,
    double dt,
    const std::vector<bool>& /*external_mask*/
) {
    const int N = static_cast<int>(system.bodies.size());
    std::vector<bool> in_subsystem(N, false);

    // --- 1. Construir el árbol de jerarquía ---
    last_root = builder.build(system);

    // --- 2. Integrar recursivamente ---
    integrate_node(*last_root, system, dt, in_subsystem);

    // --- 3. Cuerpos libres restantes → campo ---
    // (Los LEAF ya se manejan en integrate_leaf, pero por seguridad
    //  cualquier cuerpo no marcado pasa por far)
    far->step(system, dt, in_subsystem);

    ++step_counter;
}

// ============================================================================
// DESPACHADOR RECURSIVO
// ============================================================================
void HierarchicalIntegrator::integrate_node(
    HierarchyNode& node,
    NBodySystem& system,
    double dt,
    std::vector<bool>& in_subsystem
) {
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
        case HierarchyNode::Type::COMPOSITE:
            integrate_composite(node, system, dt, in_subsystem);
            break;
    }
}

// ============================================================================
// LEAF: marcar como libre (far lo integrará)
// ============================================================================
void HierarchicalIntegrator::integrate_leaf(
    HierarchyNode& node,
    NBodySystem& /*system*/,
    double /*dt*/,
    std::vector<bool>& in_subsystem
) {
    // No marcamos in_subsystem → far->step() lo integrará
    // Esto es correcto: LEAF se deja al integrador de campo
    (void)node;
    (void)in_subsystem;
    // Nada que hacer aquí: far->step() usa la máscara in_subsystem
    // y como LEAF no está marcado, lo integrará
}

// ============================================================================
// PAIR_KS: KS simple o perturbado según tidal_parameter
// ============================================================================
void HierarchicalIntegrator::integrate_pair_ks(
    HierarchyNode& node,
    NBodySystem& system,
    double dt
) {
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

    if (node.tidal_parameter < tidal_threshold) {
        // KS simple: binaria prácticamente aislada
        ks_simple.integrate(state, dt);
    } else {
        // KS perturbado: campo externo significativo
        ks_perturbed.integrate_perturbed(state, dt, system, i, j);
    }

    state.write_back(system.bodies[i], system.bodies[j]);
}

// ============================================================================
// TRIPLE_CHAIN: Chain3Integrator
// ============================================================================
void HierarchicalIntegrator::integrate_triple_chain(
    HierarchyNode& node,
    NBodySystem& system,
    double dt
) {
    const int i = node.body_indices[0];
    const int j = node.body_indices[1];
    const int k = node.body_indices[2];

    if (logger) {
        logger->log({static_cast<int>(step_counter), i, j, "ENTER_CHAIN3"});
    }

    Chain3State state = chain3.initialize(system, i, j, k);

    double t_target  = state.cm_time + dt;
    double t_achieved = 0.0;
    IntegrationParams params;
    params.abs_tol  = 1e-10;
    params.min_dtau = 1e-8;
    params.max_dtau = 1e-1;

    chain3.integrate(state, t_target, t_achieved, params, system);
    chain3.write_back(state, system, i, j, k);

    if (logger) {
        logger->log({static_cast<int>(step_counter), i, j, "EXIT_CHAIN3"});
    }
}

// ============================================================================
// COMPOSITE: integrar hijos recursivamente
// ============================================================================
void HierarchicalIntegrator::integrate_composite(
    HierarchyNode& node,
    NBodySystem& system,
    double dt,
    std::vector<bool>& in_subsystem
) {
    for (auto& child : node.children) {
        integrate_node(*child, system, dt, in_subsystem);
    }
}