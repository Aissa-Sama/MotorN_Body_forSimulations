#pragma once
#include <memory>
#include <vector>
#include <map>
#include <tuple>
#include "integrator.h"
#include "hierarchy_node.h"
#include "hierarchy_builder.h"
#include "chain3_integrator.h"
#include "archain3_integrator.h"
#include "archain3_state.h"
#include "ks_perturbed_integrator.h"
#include "ks_integrator.h"
#include "regime_logger.h"
#include "nbody_system.h"

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

    const HierarchyNode* last_tree() const { return last_root.get(); }

private:
    using TripleKey = std::tuple<int,int,int>;

    static TripleKey make_key(int i, int j, int k) {
        int a = i, b = j, c = k;
        if (a > b) std::swap(a, b);
        if (b > c) std::swap(b, c);
        if (a > b) std::swap(a, b);
        return {a, b, c};
    }

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

    void integrate_triple_ar_chain(
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

    void collect_active_ar_keys(
        const HierarchyNode& node,
        std::vector<TripleKey>& keys
    ) const;

    std::unique_ptr<Integrator>   far;
    KSIntegrator                  ks_simple;
    KSPerturbedIntegrator         ks_perturbed;
    Chain3Integrator              chain3;
    ARChain3Integrator            ar_chain_;
    HierarchyBuilder              builder;
    double                        tidal_threshold;
    RegimeLogger*                 logger;
    std::size_t                   step_counter = 0;

    std::unique_ptr<HierarchyNode> last_root;

    std::map<TripleKey, ARChain3State> ar_chain_states_;
};
