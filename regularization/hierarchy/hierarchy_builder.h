// regularization/hierarchy/hierarchy_builder.h
// FASE 6A: Añadidos PNParams dentro de Params para activación automática de PN.
#pragma once
#include <vector>
#include <tuple>
#include <map>
#include <string>
#include <memory>
#include "hierarchy_node.h"
#include "nbody_system.h"

// ============================================================================
// HIERARCHY BUILDER
//
// Construye el árbol de jerarquía a partir del estado físico del sistema.
//
// FASE 6A — ACTIVACIÓN AUTOMÁTICA DE PN:
//   Cuando PNParams::enabled=true, build() marca pn_active=true en cualquier
//   nodo cuya separación mínima cae bajo el umbral físico:
//
//     r_PN = G · m_subgrupo / (ε_PN · c²)
//
//   El HierarchicalIntegrator lee pn_active en cada nodo y selecciona
//   ARChainNPNBSIntegrator cuando corresponde.
//
//   HISTERESIS: Para evitar flip-flopping cuando sep_min ~ r_PN, el builder
//   usa un cache de histeresis (pasado desde el integrador):
//     Activar:   sep_min <  r_PN
//     Mantener:  sep_min <  2·r_PN  (si ya estaba activo)
//     Desactivar: sep_min >= 2·r_PN
// ============================================================================
class HierarchyBuilder {
public:
    // -----------------------------------------------------------------------
    // PARÁMETROS DE CONSTRUCCIÓN
    // -----------------------------------------------------------------------
    struct Params {
        double r_ks_threshold      = 1.0;
        double tidal_threshold     = 0.1;
        double strong_coupling_eta = 3.0;
        double ar_chain_threshold  = 0.5;
        double ar_chain_eta        = 1e-3;

        // ── FASE 6A — Parámetros de activación PN ───────────────────────────
        struct PNParams {
            bool   enabled             = false;  ///< false → PN nunca se activa
            double c_speed             = 1e4;    ///< Velocidad de la luz en unidades N-body
            int    pn_order            = 1;      ///< Bitmask: 1=PN1, 2=PN2, 4=PN25, 7=todos
            double activation_epsilon  = 1e-3;   ///< ε_PN: nivel de importancia relativa PN1
            double hysteresis_factor   = 2.0;    ///< Mantener PN activo hasta sep > factor×r_PN
            double bs_eps              = 1e-10;  ///< Tolerancia del GBS cuando PN activo
            double eta_pn              = 1e-3;   ///< Parámetro η del integrador PN

            /// Calcula el umbral de separación para el grupo con masa total m_total.
            /// r_PN = m_total / (ε_PN · c²)   [G=1]
            double r_pn_threshold(double m_total) const {
                if (c_speed < 1e-30 || activation_epsilon < 1e-30) return 0.0;
                return m_total / (activation_epsilon * c_speed * c_speed);
            }
        } pn;
        // ────────────────────────────────────────────────────────────────────
    };

    explicit HierarchyBuilder(const Params& p = Params{}) : params(p) {}

    // -----------------------------------------------------------------------
    // INTERFAZ PRINCIPAL
    //
    // pn_cache: mapa de histeresis. La clave es la clave canónica del subgrupo
    //   (índices ordenados, concatenados con '_'). El valor es true si el
    //   subgrupo tenía PN activo en el paso anterior.
    //   El HierarchicalIntegrator gestiona este mapa entre pasos.
    // -----------------------------------------------------------------------

    /// Construye el árbol. pn_cache es lectura/escritura: read para histeresis,
    /// write para actualizar el estado de activación de cada subgrupo.
    std::unique_ptr<HierarchyNode> build(
        const NBodySystem& system,
        std::map<std::string, bool>* pn_cache = nullptr) const;

    // -----------------------------------------------------------------------
    // MÉTODOS PÚBLICOS (visibles para tests)
    // -----------------------------------------------------------------------
    std::vector<std::pair<int,int>> find_close_pairs(const NBodySystem& system) const;
    double compute_binding_energy(const NBodySystem& system, int i, int j) const;
    double compute_tidal_parameter(const NBodySystem& system, int i, int j) const;
    std::vector<std::pair<int,int>> select_bound_pairs(const NBodySystem& system) const;
    std::vector<std::tuple<int,int,int>> find_triples(
        const std::vector<std::pair<int,int>>& pairs) const;
    bool triple_is_strongly_coupled(const NBodySystem& system, int i, int j, int k) const;
    double compute_sep_min(const NBodySystem& system, int i, int j, int k) const;

    // ── FASE 6A: helpers de clave y umbral ──────────────────────────────────

    /// Genera la clave canónica de un subgrupo a partir de sus índices ordenados.
    /// Ejemplo: índices {2, 0, 1} → "0_1_2"
    static std::string make_group_key(std::vector<int> indices);

    /// Decide si un subgrupo debe tener PN activo, considerando el cache de histeresis.
    /// Lee pn_cache[key] para el estado previo; escribe el nuevo estado.
    bool decide_pn_active(
        const std::string& key,
        double sep_min,
        double m_total,
        std::map<std::string, bool>* pn_cache) const;

private:
    Params params;
};
