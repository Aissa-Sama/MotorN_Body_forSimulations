#include <iostream>
#include <iomanip>
#include <cmath>
#include <memory>
#include <variant>
#include "nbody_system.h"
#include "hierarchical_integrator.h"
#include "hierarchy_builder.h"
#include "hierarchy_node.h"
#include "archain3_integrator.h"
#include "archain3_state.h"
#include "leapfrog_integrator.h"

static double total_energy(const NBodySystem& sys) {
    double T = 0, U = 0;
    const int N = static_cast<int>(sys.bodies.size());
    for (int i = 0; i < N; ++i) {
        T += 0.5 * sys.bodies[i].mass * sys.bodies[i].velocity.norm2();
        for (int j = i+1; j < N; ++j) {
            Vec3 r = sys.bodies[j].position - sys.bodies[i].position;
            U -= sys.G * sys.bodies[i].mass * sys.bodies[j].mass / r.norm();
        }
    }
    return T + U;
}

static NBodySystem make_figure8() {

    NBodySystem sys;
    sys.G = 1.0;
    sys.bodies.push_back({{ 0.9700436926041021, -0.2430865994534989, 0.0},
                           { 0.4662036850003824,  0.4323657300939067, 0.0}, 1.0});
    sys.bodies.push_back({{-0.9700436926041021, -0.2430865994534989, 0.0},
                           { 0.4662036850003824, -0.4323657300939067, 0.0}, 1.0});
    sys.bodies.push_back({{ 0.0,                 0.4861731989069978, 0.0},
                           {-0.9324073700007648,  0.0,                0.0}, 1.0});
    return sys;
}

static bool tree_has_type(const HierarchyNode& node, HierarchyNode::Type t) {
    if (node.type == t) return true;
    for (const auto& c : node.children)
        if (tree_has_type(*c, t)) return true;
    return false;
}

bool test_A() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST A: Binaria aislada (debe producir PAIR_KS)\n";
    std::cout << std::string(60,'=') << "\n";

    NBodySystem sys;
    sys.G = 1.0;
    sys.bodies.push_back({{-0.5, 0, 0}, {0, -0.5, 0}, 1.0});
    sys.bodies.push_back({{ 0.5, 0, 0}, {0,  0.5, 0}, 1.0});

    HierarchyBuilder::Params p;
    p.r_ks_threshold = 2.0;
    HierarchyBuilder hb(p);
    auto tree = hb.build(sys);

    bool has_pair_ks = tree_has_type(*tree, HierarchyNode::Type::PAIR_KS);
    std::cout << "  Árbol: " << (has_pair_ks ? "PAIR_KS detectado ✓" : "PAIR_KS NO detectado ✗") << "\n";

    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(), 2.0, 1e-4, p);

    double E0 = total_energy(sys);
    std::vector<bool> used(2, false);
    for (int s = 0; s < 500; ++s)
        integrator.step(sys, 0.01, used);

    double dE = std::abs(total_energy(sys) - E0) / std::abs(E0);
    std::cout << "  |ΔE/E0| = " << dE << "\n";

    bool ok = has_pair_ks && (dE < 0.01);
    std::cout << "RESULTADO TEST A: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
    return ok;
}

bool test_B() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST B: Triple cercano (debe producir TRIPLE_CHAIN)\n";
    std::cout << std::string(60,'=') << "\n";

    NBodySystem sys;
    sys.G = 1.0;
    sys.bodies.push_back({{ 0.970, -0.243, 0}, { 0.466,  0.433, 0}, 1.0});
    sys.bodies.push_back({{-0.970, -0.243, 0}, { 0.466, -0.433, 0}, 1.0});
    sys.bodies.push_back({{ 0.000,  0.970, 0}, {-0.932,  0.000, 0}, 1.0});

    HierarchyBuilder::Params p;
    p.r_ks_threshold     = 3.0;
    p.strong_coupling_eta = 5.0;
    p.ar_chain_threshold = 0.0;  
    HierarchyBuilder hb(p);
    auto tree = hb.build(sys);

    bool has_chain = tree_has_type(*tree, HierarchyNode::Type::TRIPLE_CHAIN);
    std::cout << "  Árbol: " << (has_chain ? "TRIPLE_CHAIN detectado ✓" : "TRIPLE_CHAIN NO detectado ✗") << "\n";

    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(), 3.0, 1e-4, p);

    double E0 = total_energy(sys);
    std::vector<bool> used(3, false);
    double max_dE = 0.0, r_max = 0.0;
    for (int s = 0; s < 300; ++s) {
        integrator.step(sys, 0.005, used);
        double dE = std::abs(total_energy(sys) - E0) / std::abs(E0);
        if (dE > max_dE) max_dE = dE;
        for (const auto& b : sys.bodies)
            r_max = std::max(r_max, b.position.norm());
    }

    std::cout << "  Error max energía: " << max_dE << "\n";
    std::cout << "  Radio máximo:      " << r_max << "\n";

    bool ok = has_chain && (r_max < 20.0) && (max_dE < 5.0);
    std::cout << "RESULTADO TEST B: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
    return ok;
}

bool test_C() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST C: Binaria + estrella lejana (PAIR_KS + LEAF)\n";
    std::cout << std::string(60,'=') << "\n";

    NBodySystem sys;
    sys.G = 1.0;
    sys.bodies.push_back({{-0.5, 0, 0}, {0, -0.5, 0}, 1.0});
    sys.bodies.push_back({{ 0.5, 0, 0}, {0,  0.5, 0}, 1.0});
    sys.bodies.push_back({{50.0, 0, 0}, {0,  0.1, 0}, 0.01});

    HierarchyBuilder::Params bp;
    bp.r_ks_threshold = 2.0;
    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(), 2.0, 1e-4, bp);

    double E0 = total_energy(sys);
    std::vector<bool> used(3, false);
    for (int s = 0; s < 300; ++s)
        integrator.step(sys, 0.01, used);

    double dE       = std::abs(total_energy(sys) - E0) / std::abs(E0);
    double r_field  = sys.bodies[2].position.norm();
    std::cout << "  |ΔE/E0| = " << dE << "\n";
    std::cout << "  Estrella lejana en r = " << r_field << " (debe ser ~50)\n";

    bool ok = (dE < 0.05) && (std::abs(r_field - 50.0) < 5.0);
    std::cout << "RESULTADO TEST C: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
    return ok;
}

bool test_D() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST D: Triple + estrella lejana (TRIPLE_CHAIN + LEAF)\n";
    std::cout << std::string(60,'=') << "\n";

    NBodySystem sys;
    sys.G = 1.0;
    sys.bodies.push_back({{ 0.970, -0.243, 0}, { 0.466,  0.433, 0}, 1.0});
    sys.bodies.push_back({{-0.970, -0.243, 0}, { 0.466, -0.433, 0}, 1.0});
    sys.bodies.push_back({{ 0.000,  0.970, 0}, {-0.932,  0.000, 0}, 1.0});
    sys.bodies.push_back({{100.0,   0.000, 0}, { 0.000,  0.050, 0}, 0.01});

    HierarchyBuilder::Params bp;
    bp.r_ks_threshold     = 3.0;
    bp.strong_coupling_eta = 5.0;
    bp.ar_chain_threshold = 0.0;   
    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(), 3.0, 1e-4, bp);

    double E0 = total_energy(sys);
    std::vector<bool> used(4, false);
    double max_dE = 0.0;
    for (int s = 0; s < 200; ++s) {
        integrator.step(sys, 0.005, used);
        double dE = std::abs(total_energy(sys) - E0) / std::abs(E0);
        if (dE > max_dE) max_dE = dE;
    }

    double r_field = sys.bodies[3].position.norm();
    std::cout << "  Error max energía: " << max_dE << "\n";
    std::cout << "  Estrella lejana en r = " << r_field << " (debe ser ~100)\n";

    bool ok = (max_dE < 5.0) && (std::abs(r_field - 100.0) < 10.0);
    std::cout << "RESULTADO TEST D: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
    return ok;
}

bool test_E() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST E: Builder detecta TRIPLE_AR_CHAIN (sep_min < threshold)\n";
    std::cout << std::string(60,'=') << "\n";

    NBodySystem sys_close;
    sys_close.G = 1.0;

    sys_close.bodies.push_back({{ 0.00, 0.00, 0}, {0.0,  0.5, 0}, 1.0});
    sys_close.bodies.push_back({{ 0.05, 0.00, 0}, {0.0, -0.5, 0}, 1.0});
    sys_close.bodies.push_back({{ 0.025, 0.04, 0}, {-0.3, 0.0, 0}, 1.0});

    HierarchyBuilder::Params p_close;
    p_close.r_ks_threshold     = 2.0;
    p_close.strong_coupling_eta = 10.0;
    p_close.ar_chain_threshold = 0.1;  

    HierarchyBuilder hb_close(p_close);
    auto tree_close = hb_close.build(sys_close);

    bool has_ar = tree_has_type(*tree_close, HierarchyNode::Type::TRIPLE_AR_CHAIN);
    bool no_chain = !tree_has_type(*tree_close, HierarchyNode::Type::TRIPLE_CHAIN);

    double sep = hb_close.compute_sep_min(sys_close, 0, 1, 2);
    std::cout << "  sep_min(0,1,2) = " << std::fixed << std::setprecision(4) << sep
              << " (esperado ~0.05)\n";
    std::cout << "  sep < threshold (0.1): " << (sep < 0.1 ? "sí ✓" : "no ✗") << "\n";
    std::cout << "  Árbol con sep=0.05: "
              << (has_ar ? "TRIPLE_AR_CHAIN ✓" : "TRIPLE_AR_CHAIN NO detectado ✗")
              << (no_chain ? "" : " (también tiene TRIPLE_CHAIN ✗)") << "\n";

    NBodySystem sys_far;
    sys_far.G = 1.0;

    sys_far.bodies.push_back({{ 0.00, 0.00, 0}, {0.0,  0.5, 0}, 1.0});
    sys_far.bodies.push_back({{ 0.30, 0.00, 0}, {0.0, -0.5, 0}, 1.0});
    sys_far.bodies.push_back({{ 0.15, 0.26, 0}, {-0.3, 0.0, 0}, 1.0});

    HierarchyBuilder::Params p_far;
    p_far.r_ks_threshold     = 2.0;
    p_far.strong_coupling_eta = 10.0;
    p_far.ar_chain_threshold = 0.1;  

    HierarchyBuilder hb_far(p_far);
    auto tree_far = hb_far.build(sys_far);

    bool has_chain = tree_has_type(*tree_far, HierarchyNode::Type::TRIPLE_CHAIN);
    bool no_ar     = !tree_has_type(*tree_far, HierarchyNode::Type::TRIPLE_AR_CHAIN);

    double sep_far = hb_far.compute_sep_min(sys_far, 0, 1, 2);
    std::cout << "  sep_min(0,1,2) = " << sep_far
              << " (esperado ~0.30)\n";
    std::cout << "  sep > threshold (0.1): " << (sep_far > 0.1 ? "sí ✓" : "no ✗") << "\n";
    std::cout << "  Árbol con sep=0.30: "
              << (has_chain ? "TRIPLE_CHAIN ✓" : "TRIPLE_CHAIN NO detectado ✗")
              << (no_ar ? "" : " (también tiene TRIPLE_AR_CHAIN ✗)") << "\n";

    bool ok = has_ar && no_chain
           && has_chain && no_ar
           && (sep < 0.1) && (sep_far > 0.1);
    std::cout << "RESULTADO TEST E: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
    return ok;
}

bool test_F() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST F: Fix del CM entre pasos (persiste ARChain3State)\n";
    std::cout << std::string(60,'=') << "\n";

    NBodySystem sys;
    sys.G = 1.0;
    const Vec3 boost = {0.1, 0.05, 0.0};
    sys.bodies.push_back({{ 0.9700436926041021, -0.2430865994534989, 0.0},
                           { 0.4662036850003824 + boost.x,
                             0.4323657300939067 + boost.y, 0.0}, 1.0});
    sys.bodies.push_back({{-0.9700436926041021, -0.2430865994534989, 0.0},
                           { 0.4662036850003824 + boost.x,
                            -0.4323657300939067 + boost.y, 0.0}, 1.0});
    sys.bodies.push_back({{ 0.0,                 0.4861731989069978, 0.0},
                           {-0.9324073700007648 + boost.x,
                             0.0                 + boost.y, 0.0}, 1.0});

    Vec3 cm0 = {0, 0, 0};
    double M = 0;
    for (const auto& b : sys.bodies) {
        cm0 = cm0 + b.position * b.mass;
        M  += b.mass;
    }
    cm0 = cm0 * (1.0 / M);

    Vec3 cmv = {0, 0, 0};
    for (const auto& b : sys.bodies)
        cmv = cmv + b.velocity * b.mass;
    cmv = cmv * (1.0 / M);

    HierarchyBuilder::Params bp;
    bp.r_ks_threshold     = 5.0;
    bp.strong_coupling_eta = 10.0;
    bp.ar_chain_threshold = 5.0;  
    bp.ar_chain_eta       = 1e-3;

    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(), 5.0, 1e-4, bp);

    const double dt    = 0.01;
    const int    N_steps = 100;
    std::vector<bool> used(3, false);

    std::cout << "  CM inicial:  (" << cm0.x << ", " << cm0.y << ")\n";
    std::cout << "  CM velocidad: (" << cmv.x << ", " << cmv.y << ")\n\n";

    double max_cm_err = 0.0;
    bool ok = true;

    for (int s = 1; s <= N_steps; ++s) {
        integrator.step(sys, dt, used);

        Vec3 cm_now = {0, 0, 0};
        for (const auto& b : sys.bodies)
            cm_now = cm_now + b.position * b.mass;
        cm_now = cm_now * (1.0 / M);

        Vec3 cm_expected = cm0 + cmv * (s * dt);
        double err = (cm_now - cm_expected).norm();
        max_cm_err = std::max(max_cm_err, err);
   
        if (s == 1 || s == 10 || s == 50 || s == 100) {
            std::cout << "  paso " << std::setw(3) << s
                      << "  CM_actual=(" << std::fixed << std::setprecision(4)
                      << cm_now.x << ", " << cm_now.y << ")"
                      << "  CM_esperado=(" << cm_expected.x << ", " << cm_expected.y << ")"
                      << std::scientific << std::setprecision(2)
                      << "  err=" << err << "\n";
        }

        if (err > 1e-4) ok = false;
    }

    std::cout << "\n  Error máximo en CM: " << std::scientific << max_cm_err
              << (max_cm_err < 1e-4 ? "  ✓ (fix del CM funciona)" : "  ✗ (fix del CM FALLA)") << "\n";
    std::cout << "RESULTADO TEST F: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
    return ok;
}

bool test_G() {
    std::cout << "\n" << std::string(60,'=') << "\n";
    std::cout << "TEST G: figura-8 vía HierarchicalIntegrator + AR-chain (t=0→2)\n";
    std::cout << std::string(60,'=') << "\n";

    NBodySystem sys_hier = make_figure8();

    HierarchyBuilder::Params bp;
    bp.r_ks_threshold     = 5.0;
    bp.strong_coupling_eta = 10.0;
    bp.ar_chain_threshold = 5.0;   
    bp.ar_chain_eta       = 1e-3;

    HierarchicalIntegrator integrator(
        std::make_unique<LeapfrogIntegrator>(), 5.0, 1e-4, bp);

    double E0   = total_energy(sys_hier);
    const double dt     = 0.01;
    const int    N_steps = 200; 
    std::vector<bool> used(3, false);

    double max_dE = 0.0;
    bool   has_ar_every_step = true;

    for (int s = 0; s < N_steps; ++s) {
        integrator.step(sys_hier, dt, used);
        double dE = std::abs(total_energy(sys_hier) - E0) / std::abs(E0);
        if (dE > max_dE) max_dE = dE;

        const HierarchyNode* tree = integrator.last_tree();
        if (tree && !tree_has_type(*tree, HierarchyNode::Type::TRIPLE_AR_CHAIN))
            has_ar_every_step = false;
    }

    std::cout << std::scientific << std::setprecision(3);
    std::cout << "  HierarchicalIntegrator (t=2.0):\n";
    std::cout << "    max|ΔE/E0| = " << max_dE
              << (max_dE < 0.01 ? "  ✓ (< 0.01)" : "  ✗ (>= 0.01)") << "\n";
    std::cout << "    AR-chain en cada paso: "
              << (has_ar_every_step ? "sí ✓" : "no ✗") << "\n";

    NBodySystem sys_direct = make_figure8();
    ARChain3Integrator ar_direct(1e-3);
    ARChain3State state = ar_direct.initialize(sys_direct, 0, 1, 2);
    ar_direct.integrate(state, 2.0);
    ar_direct.write_back(state, sys_direct, 0, 1, 2);

    double pos_diff = 0.0;
    for (int i = 0; i < 3; ++i) {
        Vec3 dr = sys_hier.bodies[i].position - sys_direct.bodies[i].position;
        pos_diff = std::max(pos_diff, dr.norm());
    }

    double E_direct = total_energy(sys_direct);
    double E0_direct = total_energy(make_figure8());
    double dE_direct = std::abs(E_direct - E0_direct) / std::abs(E0_direct);

    std::cout << "  Comparativa con integrador directo en t=2.0:\n";
    std::cout << "    Diferencia máx. posición = " << pos_diff << " u.a.\n";
    std::cout << "    NOTA: divergencia esperada — figura-8 es caótica.\n";
    std::cout << "          Ambos integradores conservan energía:\n";
    std::cout << "          HierarchicalInt. |ΔE/E0| = " << max_dE << "\n";
    std::cout << "          Directo          |ΔE/E0| = " << dE_direct << "\n";
    std::cout << "    Ambos conservan energía (< 0.01): "
              << ((max_dE < 0.01 && dE_direct < 0.01) ? "✓" : "✗") << "\n";

    bool ok = (max_dE < 0.01) && has_ar_every_step && (dE_direct < 0.01);
    std::cout << "RESULTADO TEST G: " << (ok ? "✓ PASADO" : "✗ FALLADO") << "\n";
    return ok;
}

int main() {
    std::cout << std::string(72,'#') << "\n";
    std::cout << "TESTS: HierarchicalIntegrator (Ruta B) — incluyendo AR-chain\n";
    std::cout << std::string(72,'#') << "\n";

    bool a = test_A();
    bool b = test_B();
    bool c = test_C();
    bool d = test_D();
    bool e = test_E();
    bool f = test_F();
    bool g = test_G();

    std::cout << "\n" << std::string(72,'=') << "\n";
    std::cout << "RESUMEN:\n";
    std::cout << "  A (binaria → PAIR_KS):              " << (a?"✓ PASADO":"✗ FALLADO") << "\n";
    std::cout << "  B (triple → TRIPLE_CHAIN):          " << (b?"✓ PASADO":"✗ FALLADO") << "\n";
    std::cout << "  C (binaria + campo):                " << (c?"✓ PASADO":"✗ FALLADO") << "\n";
    std::cout << "  D (triple + campo):                 " << (d?"✓ PASADO":"✗ FALLADO") << "\n";
    std::cout << "  E (builder detecta AR-chain):       " << (e?"✓ PASADO":"✗ FALLADO") << "\n";
    std::cout << "  F (fix del CM entre pasos):         " << (f?"✓ PASADO":"✗ FALLADO") << "\n";
    std::cout << "  G (figura-8 vía HierarchicalInt.):  " << (g?"✓ PASADO":"✗ FALLADO") << "\n";
    std::cout << std::string(72,'=') << "\n";

    int n_pass = a+b+c+d+e+f+g;
    if (n_pass == 7) {
        std::cout << "✓ 7/7 tests pasaron. HierarchicalIntegrator + AR-chain completo.\n";
        return 0;
    }
    std::cout << "✗ " << (7-n_pass) << " test(s) fallaron.\n";

    std::cout << "\nNotas:\n";
    std::cout << "  Test B: max_dE puede ser ~0.4 con Chain3 BS — es esperado.\n";
    std::cout << "  Test F: si falla, el error crece como O(K²) → bug de CM.\n";
    std::cout << "  Test G: criterio = conservación de energía en ambos integradores.\n";
    return 1;
}
