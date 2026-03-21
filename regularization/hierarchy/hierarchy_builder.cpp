// regularization/hierarchy/hierarchy_builder.cpp
// FASE 7A: Deteccion de grupos N>=4 (cuadruples, quintuples).
//          make_triple_ar_chain reemplazado por make_group_ar_chain.
#include "hierarchy_builder.h"
#include <algorithm>
#include <functional>
#include <cmath>
#include <limits>
#include <numeric>

// ============================================================================
// PASO 1: Pares cercanos
// ============================================================================
std::vector<std::pair<int,int>>
HierarchyBuilder::find_close_pairs(const NBodySystem& system) const {
    std::vector<std::pair<int,int>> pairs;
    const int N = static_cast<int>(system.bodies.size());
    for (int i = 0; i < N; ++i)
        for (int j = i+1; j < N; ++j) {
            Vec3 dr = system.bodies[j].position - system.bodies[i].position;
            if (norm(dr) < params.r_ks_threshold)
                pairs.push_back({i, j});
        }
    return pairs;
}

// ============================================================================
// PASO 2: Energía de ligadura
// ============================================================================
double HierarchyBuilder::compute_binding_energy(
    const NBodySystem& system, int i, int j) const
{
    const auto& a = system.bodies[i];
    const auto& b = system.bodies[j];
    Vec3 r = b.position - a.position;
    Vec3 v = b.velocity - a.velocity;
    double mu = (a.mass * b.mass) / (a.mass + b.mass);
    return 0.5 * mu * dot(v, v)
         - system.G * a.mass * b.mass / norm(r);
}

// ============================================================================
// PASO 3: Parámetro de marea
// ============================================================================
double HierarchyBuilder::compute_tidal_parameter(
    const NBodySystem& system, int i, int j) const
{
    const int N = static_cast<int>(system.bodies.size());
    const auto& a = system.bodies[i];
    const auto& b = system.bodies[j];
    Vec3 r_int = b.position - a.position;
    double d_int = norm(r_int);
    double F_int = system.G * a.mass * b.mass / (d_int * d_int);
    if (F_int < 1e-30) return 1e30;
    Vec3 cm = (a.mass * a.position + b.mass * b.position) / (a.mass + b.mass);
    Vec3 F_ext_vec{0, 0, 0};
    for (int k = 0; k < N; ++k) {
        if (k == i || k == j) continue;
        Vec3 r_ext = system.bodies[k].position - cm;
        double d_ext = norm(r_ext);
        double f = system.G * (a.mass + b.mass) * system.bodies[k].mass
                   / (d_ext * d_ext);
        F_ext_vec = F_ext_vec + (f / d_ext) * r_ext;
    }
    return norm(F_ext_vec) / F_int;
}

// ============================================================================
// PASO 4: Greedy de pares ligados
// ============================================================================
std::vector<std::pair<int,int>>
HierarchyBuilder::select_bound_pairs(const NBodySystem& system) const {
    auto close = find_close_pairs(system);
    struct Candidate { int i, j; double binding_energy; };
    std::vector<Candidate> candidates;
    for (auto [i, j] : close) {
        double E = compute_binding_energy(system, i, j);
        if (E < 0.0) candidates.push_back({i, j, E});
    }
    std::sort(candidates.begin(), candidates.end(),
        [](const Candidate& a, const Candidate& b){
            return a.binding_energy < b.binding_energy; });
    const int N = static_cast<int>(system.bodies.size());
    std::vector<bool> used(N, false);
    std::vector<std::pair<int,int>> selected;
    for (const auto& c : candidates)
        if (!used[c.i] && !used[c.j]) {
            used[c.i] = used[c.j] = true;
            selected.push_back({c.i, c.j});
        }
    return selected;
}

// ============================================================================
// PASO 5: Triples (legacy — para TRIPLE_CHAIN)
// ============================================================================
std::vector<std::tuple<int,int,int>>
HierarchyBuilder::find_triples(
    const std::vector<std::pair<int,int>>& pairs) const
{
    std::vector<std::tuple<int,int,int>> triples;
    const int M = static_cast<int>(pairs.size());
    for (int a = 0; a < M; ++a)
        for (int b = a+1; b < M; ++b) {
            auto [ai, aj] = pairs[a];
            auto [bi, bj] = pairs[b];
            if      (ai == bi) triples.push_back({aj, ai, bj});
            else if (ai == bj) triples.push_back({aj, ai, bi});
            else if (aj == bi) triples.push_back({ai, aj, bj});
            else if (aj == bj) triples.push_back({ai, aj, bi});
        }
    return triples;
}

bool HierarchyBuilder::triple_is_strongly_coupled(
    const NBodySystem& system, int i, int j, int k) const
{
    const auto& bi = system.bodies[i];
    const auto& bj = system.bodies[j];
    const auto& bk = system.bodies[k];
    double r_binary = norm(bj.position - bi.position);
    Vec3 cm_ij = (bi.mass * bi.position + bj.mass * bj.position)
                 / (bi.mass + bj.mass);
    double r_3cm = norm(bk.position - cm_ij);
    return r_3cm < params.strong_coupling_eta * r_binary;
}

double HierarchyBuilder::compute_sep_min(
    const NBodySystem& system, int i, int j, int k) const
{
    double r12 = norm(system.bodies[j].position - system.bodies[i].position);
    double r23 = norm(system.bodies[k].position - system.bodies[j].position);
    double r13 = norm(system.bodies[k].position - system.bodies[i].position);
    return std::min({r12, r23, r13});
}

// ============================================================================
// FASE 6A: helpers de clave y umbral PN
// ============================================================================
std::string HierarchyBuilder::make_group_key(std::vector<int> indices) {
    std::sort(indices.begin(), indices.end());
    std::string key;
    for (int i = 0; i < (int)indices.size(); ++i) {
        if (i > 0) key += '_';
        key += std::to_string(indices[i]);
    }
    return key;
}

bool HierarchyBuilder::decide_pn_active(
    const std::string& key,
    double sep_min,
    double m_total,
    std::map<std::string, bool>* pn_cache) const
{
    if (!params.pn.enabled) return false;
    const double r_on  = params.pn.r_pn_threshold(m_total);
    const double r_off = r_on * params.pn.hysteresis_factor;
    bool was_active = false;
    if (pn_cache) {
        auto it = pn_cache->find(key);
        if (it != pn_cache->end()) was_active = it->second;
    }
    bool active = was_active
                ? (sep_min < r_off)
                : (sep_min < r_on);
    if (pn_cache) (*pn_cache)[key] = active;
    return active;
}

// ============================================================================
// BUILD — construcción del árbol completo
//
// FASE 7A: Detección de grupos AR-chain de N >= 3 cuerpos.
//
// ALGORITMO:
//   1. Encontrar todos los pares dentro de r_ks_threshold.
//   2. Para N=3: buscar triples donde los 3 pares están dentro del umbral.
//   3. NUEVO (Fase 7A): Para N=4,5,...: buscar grupos donde TODOS los pares
//      del grupo están dentro del umbral y la energía colectiva < 0.
//   4. Los grupos encontrados se convierten en GROUP_AR_CHAIN.
//   5. Los pares sobrantes → PAIR_KS.
//   6. Los cuerpos sueltos → LEAF.
//
// ORDEN DE PRIORIDAD:
//   Grupos más grandes primero (más ligados). Si un cuádruple y un triple
//   comparten 3 cuerpos, el cuádruple gana (menor E_binding).
// ============================================================================
std::unique_ptr<HierarchyNode>
HierarchyBuilder::build(const NBodySystem& system,
                        std::map<std::string, bool>* pn_cache) const
{
    const int N = static_cast<int>(system.bodies.size());
    std::vector<bool> used(N, false);
    std::vector<std::unique_ptr<HierarchyNode>> nodes;

    // ── Conjunto de pares cercanos para lookup O(1) ─────────────────────────
    auto close = find_close_pairs(system);
    auto in_close = [&](int a, int b) -> bool {
        for (auto [ci, cj] : close)
            if ((ci==a && cj==b) || (ci==b && cj==a)) return true;
        return false;
    };

    // ── Helper: energía colectiva de un grupo en su CM ───────────────────────
    auto group_energy = [&](const std::vector<int>& g) -> double {
        double M = 0.0;
        Vec3 vcm;
        for (int idx : g) { M += system.bodies[idx].mass; }
        for (int idx : g)
            vcm = vcm + system.bodies[idx].velocity * (system.bodies[idx].mass / M);
        double T = 0.0, U = 0.0;
        for (int idx : g) {
            Vec3 dv = system.bodies[idx].velocity - vcm;
            T += 0.5 * system.bodies[idx].mass * dot(dv, dv);
        }
        for (int a = 0; a < (int)g.size(); ++a)
            for (int b = a+1; b < (int)g.size(); ++b) {
                Vec3 dr = system.bodies[g[b]].position - system.bodies[g[a]].position;
                U -= system.G * system.bodies[g[a]].mass * system.bodies[g[b]].mass / norm(dr);
            }
        return T + U;
    };

    // ── Helper: sep_min dentro de un grupo ───────────────────────────────────
    auto group_sep_min = [&](const std::vector<int>& g) -> double {
        double s = std::numeric_limits<double>::max();
        for (int a = 0; a < (int)g.size(); ++a)
            for (int b = a+1; b < (int)g.size(); ++b) {
                double d = norm(system.bodies[g[b]].position - system.bodies[g[a]].position);
                if (d < s) s = d;
            }
        return s;
    };

    // ── Helper: todos los pares del grupo están dentro del umbral ────────────
    auto all_pairs_close = [&](const std::vector<int>& g) -> bool {
        for (int a = 0; a < (int)g.size(); ++a)
            for (int b = a+1; b < (int)g.size(); ++b)
                if (!in_close(g[a], g[b])) return false;
        return true;
    };

    // =========================================================================
    // PASO A: Detectar grupos AR-chain de tamaño N_max hacia abajo hasta N=3
    //
    // Para cada tamaño de grupo (de mayor a menor), buscar grupos donde:
    //   - Todos los pares están dentro de r_ks_threshold
    //   - Energía colectiva < 0 (grupo ligado)
    //   - Ningún cuerpo ya asignado a otro grupo
    //
    // Cuando se encuentran múltiples grupos candidatos del mismo tamaño,
    // tomar el más ligado (menor E_binding).
    // =========================================================================

    // Máximo tamaño de grupo a detectar: N cuerpos totales (sin sentido > N)
    const int max_group = std::min(N, 6);  // limitar a 6 por costo O(N^6)

    struct GroupCandidate {
        std::vector<int> indices;
        double binding_energy;
    };

    // Iterar de grupos grandes a pequeños
    for (int gsize = max_group; gsize >= 3; --gsize) {
        // Generar todas las combinaciones de gsize cuerpos no usados
        // Algoritmo: combinaciones sin repetición con bitmask para N pequeño
        std::vector<int> free;
        for (int i = 0; i < N; ++i) if (!used[i]) free.push_back(i);
        if ((int)free.size() < gsize) continue;

        std::vector<GroupCandidate> candidates;

        // Generar combinaciones de tamaño gsize de los cuerpos libres
        std::vector<int> combo(gsize);
        std::function<void(int,int)> gen = [&](int start, int depth) {
            if (depth == gsize) {
                if (!all_pairs_close(combo)) return;
                double E = group_energy(combo);
                if (E >= 0.0) return;
                candidates.push_back({combo, E});
                return;
            }
            for (int i = start; i < (int)free.size(); ++i) {
                combo[depth] = free[i];
                gen(i+1, depth+1);
            }
        };
        gen(0, 0);

        if (candidates.empty()) continue;

        // Ordenar por energía (más ligado primero) y asignar sin solapamiento
        std::sort(candidates.begin(), candidates.end(),
            [](const GroupCandidate& a, const GroupCandidate& b){
                return a.binding_energy < b.binding_energy; });

        for (auto& cand : candidates) {
            // Verificar que ningún índice ya fue asignado en esta ronda
            bool ok = true;
            for (int idx : cand.indices) if (used[idx]) { ok = false; break; }
            if (!ok) continue;

            double sep_min = group_sep_min(cand.indices);
            double m_total = 0.0;
            for (int idx : cand.indices) m_total += system.bodies[idx].mass;

            bool pn = decide_pn_active(
                make_group_key(cand.indices), sep_min, m_total, pn_cache);

            if (sep_min < params.ar_chain_threshold) {
                // GROUP_AR_CHAIN (TTL+GBS)
                nodes.push_back(HierarchyNode::make_group_ar_chain(
                    cand.indices, cand.binding_energy, pn));
            } else {
                // Sep moderada pero grupo ligado: si es triple → TRIPLE_CHAIN
                // Si es cuádruple+ → también GROUP_AR_CHAIN (no hay Chain4)
                if (gsize == 3) {
                    nodes.push_back(HierarchyNode::make_triple_chain(
                        cand.indices[0], cand.indices[1], cand.indices[2],
                        cand.binding_energy, pn));
                } else {
                    nodes.push_back(HierarchyNode::make_group_ar_chain(
                        cand.indices, cand.binding_energy, pn));
                }
            }
            for (int idx : cand.indices) used[idx] = true;
        }
    }

    // =========================================================================
    // PASO B: Pares ligados restantes → PAIR_KS
    // =========================================================================
    struct PairCandidate { int i, j; double E; };
    std::vector<PairCandidate> all_bound;
    for (auto [i, j] : close) {
        if (used[i] || used[j]) continue;
        double E = compute_binding_energy(system, i, j);
        if (E < 0.0) all_bound.push_back({i, j, E});
    }
    std::sort(all_bound.begin(), all_bound.end(),
        [](const PairCandidate& a, const PairCandidate& b){ return a.E < b.E; });
    for (const auto& c : all_bound) {
        if (used[c.i] || used[c.j]) continue;
        double eps = compute_tidal_parameter(system, c.i, c.j);
        double m_total = system.bodies[c.i].mass + system.bodies[c.j].mass;
        double sep = norm(system.bodies[c.j].position - system.bodies[c.i].position);
        bool pn = decide_pn_active(
            make_group_key({c.i, c.j}), sep, m_total, pn_cache);
        nodes.push_back(HierarchyNode::make_pair_ks(c.i, c.j, c.E, eps, pn));
        used[c.i] = used[c.j] = true;
    }

    // =========================================================================
    // PASO C: Cuerpos libres → LEAF
    // =========================================================================
    for (int i = 0; i < N; ++i)
        if (!used[i])
            nodes.push_back(HierarchyNode::make_leaf(i));

    // =========================================================================
    // PASO D: Árbol
    // =========================================================================
    if (nodes.size() == 1) return std::move(nodes[0]);
    return HierarchyNode::make_composite(std::move(nodes));
}
