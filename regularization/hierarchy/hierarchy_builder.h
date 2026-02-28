// regularization/hierarchy/hierarchy_builder.h
#pragma once
#include <vector>
#include <tuple>
#include <memory>
#include "hierarchy_node.h"
#include "nbody_system.h"

// ============================================================================
// HIERARCHY BUILDER
//
// Construye el árbol de jerarquía a partir del estado físico del sistema.
//
// Algoritmo (según Mikkola & Aarseth, parte4_jerarquica.tex):
//   1. find_close_pairs()       → pares con r_ij < r_threshold
//   2. filter_bound()           → pares con E_bind < 0
//   3. greedy_selection()       → sin solapamientos, más ligados primero
//   4. find_triples()           → dos pares que comparten un cuerpo
//   5. classify_triple()        → Chain3 si fuertemente acoplado, KS+leaf si no
//   6. compute_tidal_parameter()→ decide KS simple vs KS perturbado
//   7. Construir árbol con los nodos resultantes
// ============================================================================
class HierarchyBuilder {
public:
    // -----------------------------------------------------------------------
    // PARÁMETROS DE CONSTRUCCIÓN
    // -----------------------------------------------------------------------
    struct Params {
        double r_ks_threshold  = 1.0;   ///< Radio para considerar pares cercanos
        double tidal_threshold = 0.1;   ///< ε < threshold → KS simple; else → KS perturbado
        double strong_coupling_eta = 3.0; ///< r_3cm < eta * r_binary → triple fuerte
    };

    explicit HierarchyBuilder(const Params& p = Params{}) : params(p) {}

    // -----------------------------------------------------------------------
    // INTERFAZ PRINCIPAL
    // -----------------------------------------------------------------------

    /// Construye el árbol completo a partir del sistema.
    std::unique_ptr<HierarchyNode> build(const NBodySystem& system) const;

    // -----------------------------------------------------------------------
    // MÉTODOS PÚBLICOS (visibles para tests)
    // -----------------------------------------------------------------------

    /// Pares con r_ij < r_ks_threshold
    std::vector<std::pair<int,int>> find_close_pairs(
        const NBodySystem& system) const;

    /// Energía de ligadura del par (i,j)
    double compute_binding_energy(
        const NBodySystem& system, int i, int j) const;

    /// Parámetro de marea: F_ext / F_int para el par (i,j)
    double compute_tidal_parameter(
        const NBodySystem& system, int i, int j) const;

    /// Pares ligados (E_bind < 0), seleccionados greedy sin solapamientos
    std::vector<std::pair<int,int>> select_bound_pairs(
        const NBodySystem& system) const;

    /// Triples: dos pares seleccionados que comparten un cuerpo
    std::vector<std::tuple<int,int,int>> find_triples(
        const std::vector<std::pair<int,int>>& pairs) const;

    /// ¿El triple (i,j,k) está fuertemente acoplado?
    bool triple_is_strongly_coupled(
        const NBodySystem& system, int i, int j, int k) const;

private:
    Params params;
};