#pragma once
#include <vector>
#include <memory>
#include <variant>
#include "vec3.h"
#include "chain_state.h"
#include "ks_state.h"
#include "archain3_state.h"  

struct HierarchyNode {

    enum class Type {
        LEAF,            
        PAIR_KS,         
        TRIPLE_CHAIN,     
        TRIPLE_AR_CHAIN,  
        COMPOSITE        
    };

    Type type = Type::LEAF;

    std::vector<int> body_indices;  

    std::vector<std::unique_ptr<HierarchyNode>> children;

    double binding_energy  = 0.0;   
    double tidal_parameter = 0.0;  

    std::variant<
        std::monostate,
        KSState,
        Chain3State,
        ARChain3State    
    > regularized_state;

    static std::unique_ptr<HierarchyNode> make_leaf(int body_idx) {
        auto n = std::make_unique<HierarchyNode>();
        n->type = Type::LEAF;
        n->body_indices = {body_idx};
        return n;
    }

    static std::unique_ptr<HierarchyNode> make_pair_ks(int i, int j,
                                                        double binding_e,
                                                        double tidal_e) {
        auto n = std::make_unique<HierarchyNode>();
        n->type = Type::PAIR_KS;
        n->body_indices = {i, j};
        n->binding_energy  = binding_e;
        n->tidal_parameter = tidal_e;
        return n;
    }

    static std::unique_ptr<HierarchyNode> make_triple_chain(int i, int j, int k,
                                                             double binding_e) {
        auto n = std::make_unique<HierarchyNode>();
        n->type = Type::TRIPLE_CHAIN;
        n->body_indices = {i, j, k};
        n->binding_energy = binding_e;
        return n;
    }

    static std::unique_ptr<HierarchyNode> make_triple_ar_chain(int i, int j, int k,
                                                                double binding_e) {
        auto n = std::make_unique<HierarchyNode>();
        n->type = Type::TRIPLE_AR_CHAIN;
        n->body_indices = {i, j, k};
        n->binding_energy = binding_e;
        return n;
    }

    static std::unique_ptr<HierarchyNode> make_composite(
        std::vector<std::unique_ptr<HierarchyNode>> child_nodes)
    {
        auto n = std::make_unique<HierarchyNode>();
        n->type = Type::COMPOSITE;
        for (const auto& c : child_nodes)
            for (int idx : c->body_indices)
                n->body_indices.push_back(idx);
        n->children = std::move(child_nodes);
        return n;
    }

    bool is_regularized() const {
        return type == Type::PAIR_KS
            || type == Type::TRIPLE_CHAIN
            || type == Type::TRIPLE_AR_CHAIN;
    }

    int size() const { return static_cast<int>(body_indices.size()); }
};
