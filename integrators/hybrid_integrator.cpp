// integrators/hybrid_integrator.cpp
// ============================================================================
// INTEGRADOR HÍBRIDO CON SOPORTE CHAIN3
//
// Jerarquía de decisión en cada paso:
//   1. ¿Hay un triple fuertemente acoplado?  → Chain3Integrator
//   2. ¿Hay pares ligados?                   → KS perturbado
//   3. El resto                              → integrador de campo (far)
//
// HOOKS para Ruta B (HierarchyBuilder):
//   detect_triple()    → HierarchyBuilder::classify()
//   handle_triple()    → HierarchyNode::TRIPLE_CHAIN::integrate()
//   handle_binary_ks() → HierarchyNode::PAIR_KS::integrate()
// ============================================================================

#include "hybrid_integrator.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include "nbody_system.h"
#include "vec3.h"

// ============================================================================
// CONSTRUCTOR
// ============================================================================
HybridIntegrator::HybridIntegrator(
    std::unique_ptr<Integrator> far_integrator,
    double r_close_,
    double ks_dt,
    RegimeLogger* logger_
)
    : far(std::move(far_integrator))
    , ks_perturbed(ks_dt)
    , chain3(1e-4)          // paso ficticio interno para Chain3
    , r_close(r_close_)
    , logger(logger_)
{}

// ============================================================================
// DETECCIÓN: ¿Tres cuerpos mutuamente cercanos?
// ============================================================================
bool HybridIntegrator::is_close_triple(
    const NBodySystem& system,
    int i, int j, int k
) const {
    const Vec3 rij = system.bodies[j].position - system.bodies[i].position;
    const Vec3 rik = system.bodies[k].position - system.bodies[i].position;
    const Vec3 rjk = system.bodies[k].position - system.bodies[j].position;
    return (norm(rij) < r_close) &&
           (norm(rik) < r_close) &&
           (norm(rjk) < r_close);
}

// ============================================================================
// DETECCIÓN: ¿Un tercer cuerpo invade la zona KS de una binaria?
// ============================================================================
bool HybridIntegrator::third_body_too_close(
    const NBodySystem& system,
    int bi, int bj,
    int& third_idx
) const {
    const int N = static_cast<int>(system.bodies.size());
    // El tercer cuerpo está "demasiado cerca" si cae dentro de r_close
    // respecto a cualquiera de los dos miembros de la binaria
    const double threshold = r_close * 3.0;  // margen más amplio para triples
    for (int k = 0; k < N; ++k) {
        if (k == bi || k == bj) continue;
        Vec3 rki = system.bodies[bi].position - system.bodies[k].position;
        Vec3 rkj = system.bodies[bj].position - system.bodies[k].position;
        if (norm(rki) < threshold || norm(rkj) < threshold) {
            third_idx = k;
            return true;
        }
    }
    return false;
}

// ============================================================================
// DETECCIÓN DE TRIPLE — HOOK para Ruta B
//
// Estrategia:
//   a) Buscar el par más cercano del sistema.
//   b) Comprobar si algún tercer cuerpo está dentro de r_close*3 de ese par.
//   c) Si sí → triple. Ordenar i,j,k para la cadena (shortest-link first).
// ============================================================================
bool HybridIntegrator::detect_triple(
    const NBodySystem& system,
    TripleCandidate& triple
) const {
    const int N = static_cast<int>(system.bodies.size());
    if (N < 3) return false;

    // Encontrar el par más cercano
    double min_dist = std::numeric_limits<double>::max();
    int best_i = -1, best_j = -1;
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double d = norm(system.bodies[j].position - system.bodies[i].position);
            if (d < min_dist) {
                min_dist = d;
                best_i = i;
                best_j = j;
            }
        }
    }
    if (min_dist >= r_close) return false;

    // ¿Hay un tercer cuerpo cerca?
    int third = -1;
    if (!third_body_too_close(system, best_i, best_j, third)) return false;

    // Ordenar los tres cuerpos para la cadena: la cadena conecta
    // los cuerpos en el orden que minimiza la suma de distancias eslabón.
    // Para 3 cuerpos hay 3 posibles órdenes: (i,j,k), (i,k,j), (j,i,k)
    // Elegimos el que minimiza max(|R1|, |R2|).
    struct Order { int a, b, c; };
    Order orders[3] = {
        {best_i, best_j, third},
        {best_i, third,  best_j},
        {best_j, best_i, third}
    };

    double best_score = std::numeric_limits<double>::max();
    Order best_order = orders[0];
    for (const auto& o : orders) {
        double d1 = norm(system.bodies[o.b].position - system.bodies[o.a].position);
        double d2 = norm(system.bodies[o.c].position - system.bodies[o.b].position);
        double score = std::max(d1, d2);
        if (score < best_score) {
            best_score = score;
            best_order = o;
        }
    }

    // Calcular energía de ligadura del triple (estimación rápida)
    const auto& a = system.bodies[best_order.a];
    const auto& b = system.bodies[best_order.b];
    const auto& c = system.bodies[best_order.c];
    Vec3 vcm = (a.mass * a.velocity + b.mass * b.velocity + c.mass * c.velocity)
               / (a.mass + b.mass + c.mass);
    double T = 0.5 * a.mass * dot(a.velocity - vcm, a.velocity - vcm)
             + 0.5 * b.mass * dot(b.velocity - vcm, b.velocity - vcm)
             + 0.5 * c.mass * dot(c.velocity - vcm, c.velocity - vcm);
    double U = -a.mass * b.mass / norm(b.position - a.position)
               -b.mass * c.mass / norm(c.position - b.position)
               -a.mass * c.mass / norm(c.position - a.position);

    triple = { best_order.a, best_order.b, best_order.c, T + U };
    return true;
}

// ============================================================================
// INTEGRACIÓN DE TRIPLE — HOOK para Ruta B
// ============================================================================
void HybridIntegrator::handle_triple(
    NBodySystem& system,
    const TripleCandidate& triple,
    double dt
) {
    if (logger) {
        logger->log({
            static_cast<int>(step_counter),
            triple.i, triple.j,
            "ENTER_CHAIN3"
        });
    }

    Chain3State state = chain3.initialize(system, triple.i, triple.j, triple.k);

    double t_target  = state.cm_time + dt;
    double t_achieved = 0.0;
    IntegrationParams params;
    params.abs_tol = 1e-10;
    params.rel_tol = 1e-10;
    params.min_dtau = 1e-8;
    params.max_dtau = 1e-1;

    chain3.integrate(state, t_target, t_achieved, params, system);
    chain3.write_back(state, system, triple.i, triple.j, triple.k);

    if (logger) {
        logger->log({
            static_cast<int>(step_counter),
            triple.i, triple.j,
            "EXIT_CHAIN3"
        });
    }
}

// ============================================================================
// INTEGRACIÓN DE BINARIA KS — HOOK para Ruta B
// ============================================================================
void HybridIntegrator::handle_binary_ks(
    NBodySystem& system,
    const BinaryPair& bin,
    double dt
) {
    if (logger) {
        logger->log({
            static_cast<int>(step_counter),
            bin.i, bin.j,
            "ENTER_KS_PERTURBED"
        });
    }

    BinaryState state(system.bodies[bin.i], system.bodies[bin.j]);
    ks_perturbed.integrate_perturbed(state, dt, system, bin.i, bin.j);
    state.write_back(system.bodies[bin.i], system.bodies[bin.j]);

    if (logger) {
        logger->log({
            static_cast<int>(step_counter),
            bin.i, bin.j,
            "EXIT_KS_PERTURBED"
        });
    }
}

// ============================================================================
// DETECCIÓN DE PARES LIGADOS (sin cambios respecto al código original)
// ============================================================================
bool HybridIntegrator::is_bound_binary(
    const NBodySystem& system,
    int i, int j,
    double& binding_energy
) const {
    const auto& a = system.bodies[i];
    const auto& b = system.bodies[j];
    Vec3 r = b.position - a.position;
    Vec3 v = b.velocity - a.velocity;
    double rnorm = norm(r);
    if (rnorm > r_close) return false;
    double mu = (a.mass * b.mass) / (a.mass + b.mass);
    double kinetic   = 0.5 * mu * dot(v, v);
    double potential = -(a.mass * b.mass) / rnorm;
    binding_energy   = kinetic + potential;
    return binding_energy < 0.0;
}

std::vector<BinaryPair>
HybridIntegrator::detect_all_bound_pairs(
    const NBodySystem& system
) const {
    std::vector<BinaryPair> all_pairs;
    const int N = static_cast<int>(system.bodies.size());
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            double E;
            if (is_bound_binary(system, i, j, E))
                all_pairs.push_back({i, j, E});
        }
    }
    return all_pairs;
}

std::vector<BinaryPair>
HybridIntegrator::select_optimal_binaries(
    std::vector<BinaryPair>& candidates,
    std::vector<bool>& used,
    const NBodySystem& system
) const {
    std::sort(candidates.begin(), candidates.end(),
        [](const BinaryPair& a, const BinaryPair& b) {
            return a.binding_energy < b.binding_energy;
        });

    double energy_threshold = (system.bodies.size() > 2) ? -1.0 : -0.01;

    std::vector<BinaryPair> selected;
    for (const auto& pair : candidates) {
        if (!used[pair.i] && !used[pair.j]) {
            if (pair.binding_energy < energy_threshold) {
                used[pair.i] = true;
                used[pair.j] = true;
                selected.push_back(pair);
            }
        }
    }
    return selected;
}

// ============================================================================
// SEPARACIÓN DE FUERZAS (sin cambios)
// ============================================================================
void HybridIntegrator::compute_external_forces(
    const NBodySystem& system,
    const std::vector<bool>& in_binary,
    std::vector<Vec3>& acc_binary,
    std::vector<Vec3>& acc_field
) const {
    const int N = static_cast<int>(system.bodies.size());
    acc_binary.assign(N, {0, 0, 0});
    acc_field.assign(N,  {0, 0, 0});

    for (int i = 0; i < N; ++i) {
        if (!in_binary[i]) continue;
        for (int j = 0; j < N; ++j) {
            if (in_binary[j]) continue;
            Vec3 r = system.bodies[j].position - system.bodies[i].position;
            double d = norm(r) + 1e-9;
            acc_binary[i] = acc_binary[i] + (system.G * system.bodies[j].mass / (d*d*d)) * r;
        }
    }
    for (int i = 0; i < N; ++i) {
        if (in_binary[i]) continue;
        for (int j = 0; j < N; ++j) {
            if (!in_binary[j]) continue;
            Vec3 r = system.bodies[j].position - system.bodies[i].position;
            double d = norm(r) + 1e-9;
            acc_field[i] = acc_field[i] + (system.G * system.bodies[j].mass / (d*d*d)) * r;
        }
    }
}

// ============================================================================
// PASO PRINCIPAL
// ============================================================================
void HybridIntegrator::step(
    NBodySystem& system,
    double dt,
    const std::vector<bool>& /*unused*/
) {
    const int N = static_cast<int>(system.bodies.size());
    std::vector<bool> in_subsystem(N, false);

    // -----------------------------------------------------------------------
    // 1. ¿HAY UN TRIPLE FUERTEMENTE ACOPLADO?
    //    Si sí: Chain3 se encarga de los tres cuerpos.
    //    HOOK Ruta B: sustituir detect_triple + handle_triple por árbol jerárquico.
    // -----------------------------------------------------------------------
    TripleCandidate triple;
    if (detect_triple(system, triple)) {
        in_subsystem[triple.i] = true;
        in_subsystem[triple.j] = true;
        in_subsystem[triple.k] = true;
        handle_triple(system, triple, dt);
    }
    else {
        // -------------------------------------------------------------------
        // 2. DETECCIÓN Y SELECCIÓN DE BINARIAS KS
        //    HOOK Ruta B: sustituir por HierarchyBuilder::select_pairs()
        // -------------------------------------------------------------------
        auto all_pairs = detect_all_bound_pairs(system);
        auto selected  = select_optimal_binaries(all_pairs, in_subsystem, system);

        // -------------------------------------------------------------------
        // 3. INTEGRAR BINARIAS CON KS PERTURBADO
        //    HOOK Ruta B: sustituir por HierarchyNode::PAIR_KS::integrate()
        // -------------------------------------------------------------------
        for (const auto& bin : selected)
            handle_binary_ks(system, bin, dt);

        // -------------------------------------------------------------------
        // 4. INTEGRAR CUERPOS LIBRES CON EL INTEGRADOR DE CAMPO
        //    in_subsystem actúa como máscara: true = ya integrado, skip
        //    HOOK Ruta B: los nodos LEAF también pasan por aquí
        // -------------------------------------------------------------------
        far->step(system, dt, in_subsystem);
    }

    ++step_counter;
}
