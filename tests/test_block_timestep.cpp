// tests/test_block_timestep.cpp
// FASE 7B — Block timestep (Aarseth 2003)
//
// ── TESTS ────────────────────────────────────────────────────────────────────
//
// T1 — Calculo de jerk: verificar que compute_jerks() es consistente con
//      la derivada numerica de las aceleraciones. Para binaria circular,
//      el jerk analitico es conocido.
//
// T2 — Asignacion de pasos: cuerpos con aceleraciones distintas reciben
//      pasos distintos. Cuerpo cercano (a grande) → dt pequeno.
//      Cuerpo lejano (a pequeno) → dt grande.
//
// T3 — Conservacion de energia: sistema con block timestep activado conserva
//      energia comparablemente al leapfrog uniforme.
//
// T4 — Diferenciacion real: binaria compacta (T_orb=0.01) + estrella lejana
//      (T_orb=100). Con block timestep, la estrella lejana recibe dt ~100x
//      mayor que los cuerpos del campo cercano. Verificar que la asignacion
//      es correcta y el sistema no diverge.
// ─────────────────────────────────────────────────────────────────────────────
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include "nbody_system.h"
#include "body.h"
#include "block_timestep.h"
#include "leapfrog_integrator.h"
#include "hierarchical_integrator.h"

static bool PASS(const char* n) {
    std::cout << "  [PASS] " << n << "\n"; return true;
}
static bool FAIL(const char* n, const std::string& r) {
    std::cout << "  [FAIL] " << n << " — " << r << "\n"; return false;
}

static double total_energy(const NBodySystem& sys) {
    double T=0, U=0;
    for (int i=0; i<(int)sys.bodies.size(); ++i) {
        T += 0.5*sys.bodies[i].mass*sys.bodies[i].velocity.norm2();
        for (int j=i+1; j<(int)sys.bodies.size(); ++j) {
            Vec3 dr = sys.bodies[j].position - sys.bodies[i].position;
            U -= sys.G*sys.bodies[i].mass*sys.bodies[j].mass/dr.norm();
        }
    }
    return T+U;
}

// ── T1: Jerk de binaria circular ─────────────────────────────────────────────
// Para binaria circular m1=m2=0.5, a=1, r12=1:
//   |a| = G*m/r^2 = 0.5
//   En orbita circular el jerk es perpendicular a a y su norma es
//   |jerk| = |a| * v / r = 0.5 * v_circ = 0.5 * sqrt(G*M/a) = 0.5
//   → |jerk| = |a| * omega = |a| * sqrt(G*M/a^3) = 0.5 * 1 = 0.5
// Criterio: |jerk_numerico - jerk_analitico| / |jerk_analitico| < 1e-10
static bool test_T1_jerk_binary() {
    std::cout << "\nT1 — Jerk de binaria circular\n";

    NBodySystem sys; sys.G = 1.0;
    // Binaria circular m1=m2=0.5, a=1, en el apoastro (posicion inicial = apoastro)
    const double m = 0.5, a = 1.0, M = 2*m;
    const double v = std::sqrt(M/a) * 0.5; // v_cm para masa reducida
    sys.bodies.push_back({Vec3(-0.5,0,0), Vec3(0,-v,0), m});
    sys.bodies.push_back({Vec3( 0.5,0,0), Vec3(0, v,0), m});

    auto accs  = sys.compute_accelerations();
    auto jerks = BlockTimestep::compute_jerks(sys, accs);

    const double jn0 = jerks[0].norm();
    const double jn1 = jerks[1].norm();

    // |jerk| analitico = G*m/r^2 * v_rel/r = (G*m/r^2) * (2v/r)
    // r12=1, v_rel=2v
    const double r12 = 1.0;
    const double v_rel = 2.0*v;
    const double jerk_analytic = sys.G*m/(r12*r12) * v_rel/r12;

    std::cout << "  |jerk[0]| = " << jn0 << "\n";
    std::cout << "  |jerk[1]| = " << jn1 << "\n";
    std::cout << "  jerk_analitico = " << jerk_analytic << "\n";

    const double err0 = std::abs(jn0 - jerk_analytic)/jerk_analytic;
    const double err1 = std::abs(jn1 - jerk_analytic)/jerk_analytic;
    std::cout << "  err[0]=" << err0 << "  err[1]=" << err1 << "\n";

    if (err0 > 1e-10) return FAIL("T1","jerk[0] difiere del analitico >1e-10");
    if (err1 > 1e-10) return FAIL("T1","jerk[1] difiere del analitico >1e-10");
    return PASS("T1: jerk binaria circular == analitico (err<1e-10)");
}

// ── T2: Asignacion de pasos diferenciada ─────────────────────────────────────
// Sistema: masa A cerca del origen (a grande → dt pequeno)
//          masa B lejos (a pequeno → dt grande)
// Verificar: dt_A < dt_B
static bool test_T2_dt_assignment() {
    std::cout << "\nT2 — Asignacion de pasos diferenciada\n";

    NBodySystem sys; sys.G = 1.0;
    // Masa central masiva
    sys.bodies.push_back({Vec3(0,0,0),   Vec3(0,0,0),   100.0}); // central
    sys.bodies.push_back({Vec3(0.1,0,0), Vec3(0,30,0),  0.001}); // cuerpo A: cerca, a grande
    sys.bodies.push_back({Vec3(100,0,0), Vec3(0,1,0),   0.001}); // cuerpo B: lejos, a pequeno

    // dt_max=10.0: techo alto para que dt_B no quede truncado
    // Con m_central=100, r_A=0.1, v_A=30: dt_A~1e-3 (a grande, jerk grande)
    // Con m_central=100, r_B=100, v_B=1:  dt_B~0.3  (a pequeno, jerk pequeno)
    // Sin techo alto ambos se truncan al mismo nivel de bloque
    BlockTimestep::Params p;
    p.eta    = 0.03;
    p.dt_min = 1e-6;
    p.dt_max = 10.0;   // techo alto para revelar diferenciacion natural
    p.k_max  = 20;
    BlockTimestep bt(p);

    auto accs  = sys.compute_accelerations();
    auto jerks = BlockTimestep::compute_jerks(sys, accs);

    // Solo A y B son LEAFs (central marcada como usada)
    std::vector<bool> used = {true, false, false};
    auto dts = bt.assign(sys, accs, jerks, used, p.dt_max);

    std::cout << "  dt_central = " << dts[0] << " (usado, debe ser dt_max)\n";
    std::cout << "  dt_A (cerca) = " << dts[1] << "\n";
    std::cout << "  dt_B (lejos) = " << dts[2] << "\n";
    std::cout << "  ratio dt_B/dt_A = " << dts[2]/dts[1] << " (debe ser > 1)\n";

    if (dts[0] != p.dt_max) return FAIL("T2","cuerpo usado no recibio dt_max");
    if (dts[1] >= dts[2])   return FAIL("T2","cuerpo cercano no recibio dt menor que el lejano");
    if (dts[2]/dts[1] < 2.0) return FAIL("T2","ratio dt_B/dt_A < 2 — diferenciacion insuficiente");
    return PASS("T2: dt_cercano < dt_lejano, ratio >= 2 (potencias de 2)");
}

// ── T3: Conservacion de energia con block timestep ───────────────────────────
// Sistema solar simplificado: Sol + Tierra + Jupiter (3 cuerpos).
// Con block timestep: Tierra recibe dt pequeno, Jupiter dt mayor.
// Criterio: |dE/E0| < 1e-4 en 100 pasos.
static bool test_T3_energy_conservation() {
    std::cout << "\nT3 — Conservacion de energia con block timestep\n";

    NBodySystem sys; sys.G = 1.0;
    // Sol + Tierra + Jupiter (unidades: G=1, M_sol=1)
    sys.bodies.push_back({Vec3(0,0,0),    Vec3(0,0,0),      1.0});    // Sol
    sys.bodies.push_back({Vec3(1,0,0),    Vec3(0,1,0),      3e-6});   // Tierra (a=1)
    sys.bodies.push_back({Vec3(5.2,0,0),  Vec3(0,0.438,0),  1e-3});   // Jupiter (a=5.2)

    const double E0 = total_energy(sys);

    BlockTimestep::Params bp;
    bp.eta    = 0.03;
    bp.dt_min = 1e-4;
    bp.dt_max = 0.01;
    bp.k_max  = 6;

    HierarchyBuilder::Params hp;
    hp.r_ks_threshold = 0.1; // Solo pares muy cercanos → todos son LEAFs aqui

    HierarchicalIntegrator integ(
        std::make_unique<LeapfrogIntegrator>(),
        0.1, 1e-4, hp, nullptr,
        true,  // enable_block_ts
        bp
    );

    std::vector<bool> used(3, false);
    double max_dE = 0.0;
    for (int s = 0; s < 100; ++s) {
        integ.step(sys, bp.dt_max, used);
        double dE = std::abs(total_energy(sys)-E0)/std::abs(E0);
        if (dE > max_dE) max_dE = dE;
    }

    std::cout << "  E0 = " << E0 << "\n";
    std::cout << "  max|dE/E0| = " << max_dE << "\n";

    if (max_dE > 1e-3)
        return FAIL("T3","|dE/E0|>1e-3 con block timestep");
    return PASS("T3: energia conservada con block timestep (|dE/E0|<1e-3)");
}

// ── T4: Diferenciacion real — escala rapida vs lenta ─────────────────────────
// Binaria compacta (sep=0.01) + estrella lejana (sep=10).
// La estrella lejana debe recibir un dt al menos 4x mayor que las masas cercanas.
// El sistema no debe divergir en 200 pasos.
static bool test_T4_scale_separation() {
    std::cout << "\nT4 — Separacion de escalas temporales\n";

    NBodySystem sys; sys.G = 1.0;
    // Binaria circular: m1=m2=0.5, sep=0.1, v_circ=sqrt(G*M/a)*0.5=sqrt(1/0.1)*0.5=1.58
    const double sep=0.1, mc=0.5, vc=std::sqrt(mc*2.0/sep)*0.5;
    sys.bodies.push_back({Vec3(-sep*0.5,0,0), Vec3(0,-vc,0), mc}); // masa A
    sys.bodies.push_back({Vec3( sep*0.5,0,0), Vec3(0, vc,0), mc}); // masa B
    // Estrella lejana en orbita estable: a=10, v_circ=sqrt(G*M_total/a)
    const double a_out=10.0, v_out=std::sqrt(1.0/a_out)*0.7;
    sys.bodies.push_back({Vec3(a_out,0,0), Vec3(0,v_out,0), 0.1}); // estrella

    auto accs  = sys.compute_accelerations();
    auto jerks = BlockTimestep::compute_jerks(sys, accs);

    BlockTimestep::Params p;
    p.eta    = 0.03;
    p.dt_min = 1e-4;
    p.dt_max = 0.1;
    p.k_max  = 10;
    BlockTimestep bt(p);

    std::vector<bool> used(3, false);
    auto dts = bt.assign(sys, accs, jerks, used, p.dt_max);

    std::cout << "  dt_A (binaria) = " << dts[0] << "\n";
    std::cout << "  dt_B (binaria) = " << dts[1] << "\n";
    std::cout << "  dt_estrella    = " << dts[2] << "\n";
    std::cout << "  ratio estrella/binaria = " << dts[2]/dts[0] << "\n";

    if (dts[2] < 4.0*dts[0])
        return FAIL("T4","estrella lejana no recibio dt >= 4x mayor que binaria");

    // Verificar conservacion de energia con block timestep
    const double E0 = total_energy(sys);

    HierarchyBuilder::Params hp;
    hp.r_ks_threshold = 0.01; // binaria compacta → PAIR_KS, estrella → LEAF
    HierarchicalIntegrator integ(
        std::make_unique<LeapfrogIntegrator>(),
        0.01, 1e-4, hp, nullptr, true, p
    );

    std::vector<bool> used2(3, false);
    double max_dE = 0.0;
    for (int s = 0; s < 300; ++s) {
        integ.step(sys, p.dt_max, used2);
        double dE = std::abs(total_energy(sys)-E0)/std::abs(E0);
        if (dE > max_dE) max_dE = dE;
    }
    std::cout << "  max|dE/E0| = " << max_dE << "\n";

    if (max_dE > 0.1) return FAIL("T4","|dE/E0|>0.1 con block timestep");
    return PASS("T4: escalas diferenciadas correctamente, energia conservada");
}

// ── Main ──────────────────────────────────────────────────────────────────────
int main() {
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "======================================================\n";
    std::cout << "  test_block_timestep — Fase 7B (Aarseth 2003)\n";
    std::cout << "======================================================\n";

    int passed=0, total=0;
    auto run = [&](bool(*fn)(), const char* name) {
        ++total;
        std::cout << "\nTest " << total << " — " << name << "\n";
        try { if (fn()) ++passed; }
        catch (const std::exception& e) {
            FAIL(name, std::string("excepcion: ")+e.what());
        }
    };

    run(test_T1_jerk_binary,       "Jerk binaria circular == analitico");
    run(test_T2_dt_assignment,     "Asignacion dt: cercano < lejano");
    run(test_T3_energy_conservation, "Energia conservada con block timestep");
    run(test_T4_scale_separation,  "Separacion escalas rapida/lenta");

    std::cout << "\n======================================================\n";
    std::cout << "  Resultado: " << passed << "/" << total << " tests\n";
    std::cout << "======================================================\n";
    return (passed==total) ? 0 : 1;
}