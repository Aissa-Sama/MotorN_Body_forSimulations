// integrators/hierarchical_integrator.h
#pragma once
#include <memory>
#include <vector>
#include "integrator.h"
#include "hierarchy_node.h"
#include "hierarchy_builder.h"
#include "chain3_integrator.h"
#include "ks_perturbed_integrator.h"
#include "ks_integrator.h"
#include "regime_logger.h"
#include "nbody_system.h"

// ============================================================================
// HIERARCHICAL INTEGRATOR — Ruta B
//
// Reemplaza HybridIntegrator usando el árbol de HierarchyBuilder.
//
// Cada paso:
//   1. HierarchyBuilder::build() → árbol del instante actual
//   2. Recorrer el árbol e integrar cada nodo con su método especializado:
//        LEAF         → far (Leapfrog/RK4)
//        PAIR_KS      → KS simple o KS perturbado (según tidal_parameter)
//        TRIPLE_CHAIN → Chain3Integrator
//        COMPOSITE    → integrar hijos recursivamente
//
// Ventajas sobre HybridIntegrator:
//   - Detección automática de triples
//   - Clasificación automática KS simple vs perturbado
//   - Estructura extensible: añadir COMPOSITE anidado en el futuro es trivial
// ============================================================================
class HierarchicalIntegrator : public Integrator {
public:
    HierarchicalIntegrator(
        std::unique_ptr<Integrator> far_integrator,
        double r_ks_threshold,
        double ks_internal_dt,
        const HierarchyBuilder::Params& builder_params = HierarchyBuilder::Params{},
        RegimeLogger* logger = nullptr
    );

    void step(
        NBodySystem& system,
        double dt,
        const std::vector<bool>& used
    ) override;

    /// Acceso al árbol del último paso (útil para tests y diagnóstico)
    const HierarchyNode* last_tree() const { return last_root.get(); }

private:
    // -----------------------------------------------------------------------
    // INTEGRACIÓN POR TIPO DE NODO
    // -----------------------------------------------------------------------

    void integrate_node(
        HierarchyNode& node,
        NBodySystem& system,
        double dt,
        std::vector<bool>& in_subsystem
    );

    void integrate_leaf(
        HierarchyNode& node,
        NBodySystem& system,
        double dt,
        std::vector<bool>& in_subsystem
    );

    void integrate_pair_ks(
        HierarchyNode& node,
        NBodySystem& system,
        double dt
    );

    void integrate_triple_chain(
        HierarchyNode& node,
        NBodySystem& system,
        double dt
    );

    void integrate_composite(
        HierarchyNode& node,
        NBodySystem& system,
        double dt,
        std::vector<bool>& in_subsystem
    );

    // -----------------------------------------------------------------------
    // MIEMBROS
    // -----------------------------------------------------------------------
    std::unique_ptr<Integrator>   far;
    KSIntegrator                  ks_simple;
    KSPerturbedIntegrator         ks_perturbed;
    Chain3Integrator              chain3;
    HierarchyBuilder              builder;
    double                        tidal_threshold;
    RegimeLogger*                 logger;
    std::size_t                   step_counter = 0;

    std::unique_ptr<HierarchyNode> last_root;  ///< Árbol del último paso
};