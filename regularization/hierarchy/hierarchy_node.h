// regularization/hierarchy/hierarchy_node.h
#pragma once
#include <vector>
#include <memory>
#include <variant>
#include "vec3.h"
#include "chain_state.h"
#include "ks_state.h"

// ============================================================================
// NODO DEL ÁRBOL DE JERARQUÍA
//
// El árbol representa la estructura dinámica del sistema en cada paso:
//
//   COMPOSITE
//   ├── TRIPLE_CHAIN  (cuerpos 0,1,2 → Chain3)
//   ├── PAIR_KS       (cuerpos 3,4   → KS perturbado)
//   └── LEAF          (cuerpo 5      → campo Leapfrog)
//
// Ruta B reemplaza los métodos detect_triple/handle_triple del
// HybridIntegrator con la construcción y recorrido de este árbol.
// ============================================================================

struct HierarchyNode {

    // -----------------------------------------------------------------------
    // TIPO DE NODO
    // -----------------------------------------------------------------------
    enum class Type {
        LEAF,           ///< Cuerpo individual → integrador de campo
        PAIR_KS,        ///< Par ligado → KS simple o perturbado
        TRIPLE_CHAIN,   ///< Triple fuertemente acoplado → Chain3
        COMPOSITE       ///< Contiene sub-nodos (raíz o subsistemas jerárquicos)
    };

    Type type = Type::LEAF;

    // -----------------------------------------------------------------------
    // ÍNDICES EN NBodySystem
    // -----------------------------------------------------------------------
    std::vector<int> body_indices;  ///< Índices de los cuerpos de este nodo

    // -----------------------------------------------------------------------
    // HIJOS (solo para COMPOSITE)
    // -----------------------------------------------------------------------
    std::vector<std::unique_ptr<HierarchyNode>> children;

    // -----------------------------------------------------------------------
    // PROPIEDADES DE ACOPLAMIENTO
    // -----------------------------------------------------------------------
    double binding_energy  = 0.0;   ///< < 0 si el subsistema está ligado
    double tidal_parameter = 0.0;   ///< F_ext / F_int — controla KS simple vs perturbado

    // -----------------------------------------------------------------------
    // ESTADO REGULARIZADO (KS o Chain3, si aplica)
    //
    // LEAF       → monostate (sin estado regularizado)
    // PAIR_KS    → KSState   (estado KS del par)
    // TRIPLE_CHAIN → Chain3State (estado chain del triple)
    // -----------------------------------------------------------------------
    std::variant<
        std::monostate,
        KSState,
        Chain3State
    > regularized_state;

    // -----------------------------------------------------------------------
    // CONSTRUCTORES DE CONVENIENCIA
    // -----------------------------------------------------------------------

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

    static std::unique_ptr<HierarchyNode> make_composite(
        std::vector<std::unique_ptr<HierarchyNode>> child_nodes)
    {
        auto n = std::make_unique<HierarchyNode>();
        n->type = Type::COMPOSITE;
        // Recopilar todos los índices de los hijos
        for (const auto& c : child_nodes) {
            for (int idx : c->body_indices)
                n->body_indices.push_back(idx);
        }
        n->children = std::move(child_nodes);
        return n;
    }

    // -----------------------------------------------------------------------
    // UTILIDADES
    // -----------------------------------------------------------------------

    bool is_regularized() const {
        return type == Type::PAIR_KS || type == Type::TRIPLE_CHAIN;
    }

    int size() const { return static_cast<int>(body_indices.size()); }
};