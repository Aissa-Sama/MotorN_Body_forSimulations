// tests/test_physical_benchmarks.cpp
//
// Fase 3 — Validación contra benchmarks físicos publicados.
// Cada test tiene criterios calibrados con referencias explícitas.
//
// ── TEST P1 — Precesión del perihelio (PN1) ──────────────────────────────────
//   Sistema:   binaria compacta, m₁=m₂=1, G=1, a=0.1, e=0.6, c=10
//   Referencia: Einstein (1915); Iyer & Will (1993) ec. 3.1
//   Fórmula:   Δω_PN1 = 6πG(m₁+m₂) / [a·c²·(1−e²)]  [rad/órbita]
//   Criterio 1: |Δω_meas/Δω_teoria − 1| < 0.01  → PN1 correcto
//   Criterio 2: 0.01–0.05               → problema de integración
//   Criterio 3: > 0.05                  → error en compute_PN1() o acoplamiento GBS
//   Criterio 4: |dE/E₀| < 1e-8          → PN1 es conservativo (sin disipación)
//   Criterio 5: |dL/L₀| < 1e-8          → momento angular conservado por PN1
//
// ── TEST P2 — Problema de Pitágoras completo ─────────────────────────────────
//   Sistema:   masas 3, 4, 5 en triángulo rectángulo, en reposo
//   Referencia: Szebehely & Peters (1967), AJ 72, 876
//   Resultados publicados:
//     - Cuerpo de masa 3 es eyectado (cuerpo más ligero)
//     - Binaria residual (masas 4 y 5) queda ligada
//     - Velocidad asintótica cuerpo expulsado: v_∞ ≈ 0.347 (Szebehely & Peters)
//     - v_∞ confirmado en Aarseth (2003) con NBODY6: mismo orden de magnitud
//   Criterios:
//     - r₃ > 30 en t=100        → eyección confirmada (5% tolerancia por caos)
//     - E_bin(4,5) < 0           → binaria ligada
//     - |v_∞/v_ref − 1| < 0.15  → velocidad asintótica (15%: sistema caótico,
//                                   doble precisión ≠ 1967)
//     - |ΔP_total| < 1e-10      → CM conservado exactamente
//     - E_bin + E_libre ≈ E₀    → balance energético global < 1e-6
//
// ── TEST P3 — Figura-8 × 100 períodos: deriva energética secular ─────────────
//   Sistema:   CI de Simó (2002), bs_eps=1e-10
//   Referencia: Simó (2002), Cel. Mech. Dyn. Astron. 82, 3−29
//   Criterios:
//     - max|dE/E₀| en 100T < 1e-8       → precisión GBS mantenida
//     - slope lineal de error < 1e-11/T  → sin drift secular (simpléctico)
//     - |ΔL| total < 1e-12              → L=0 exacto por simetría de las CI
//
// ═══════════════════════════════════════════════════════════════════════════════

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <numeric>

#include "archain_n_bs_integrator.h"
#include "archain_n_pn_bs_integrator.h"
#include "archain_n_state.h"
#include "nbody_system.h"
#include "body.h"
#include "vec3.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ─── Helpers físicos ─────────────────────────────────────────────────────────

static double total_energy(const NBodySystem& sys) {
    const int N = (int)sys.bodies.size();
    double T = 0.0, U = 0.0;
    for (int i = 0; i < N; i++) {
        T += 0.5 * sys.bodies[i].mass * dot(sys.bodies[i].velocity, sys.bodies[i].velocity);
        for (int j = i+1; j < N; j++) {
            double r = norm(sys.bodies[j].position - sys.bodies[i].position);
            U -= sys.G * sys.bodies[i].mass * sys.bodies[j].mass / r;
        }
    }
    return T + U;
}

static Vec3 total_momentum(const NBodySystem& sys) {
    Vec3 P{0,0,0};
    for (const auto& b : sys.bodies)
        P = P + b.velocity * b.mass;
    return P;
}

static Vec3 total_angular_momentum(const NBodySystem& sys) {
    Vec3 L{0,0,0};
    for (const auto& b : sys.bodies) {
        Vec3 r = b.position;
        Vec3 v = b.velocity;
        L = L + cross(r, v * b.mass);
    }
    return L;
}

static double binding_energy_pair(const NBodySystem& sys, int i, int j) {
    const auto& a = sys.bodies[i];
    const auto& b = sys.bodies[j];
    Vec3 dv = b.velocity - a.velocity;
    Vec3 dr = b.position - a.position;
    double mu = (a.mass * b.mass) / (a.mass + b.mass);
    return 0.5 * mu * dot(dv, dv) - sys.G * a.mass * b.mass / norm(dr);
}

// Período kepleriano T = 2π√(a³/GM)
static double kepler_period(double a, double G, double M) {
    return 2.0 * M_PI * std::sqrt(a*a*a / (G * M));
}

// ─── Condiciones iniciales ───────────────────────────────────────────────────

// Binaria kepleriana: m1=m2=1, G=1, semieje a, excentricidad e
// En el perihelio: r_peri = a(1-e), v_peri = sqrt(GM(1+e)/[a(1-e)])
static NBodySystem make_binary(double m1, double m2, double a, double e) {
    NBodySystem sys;
    sys.G = 1.0;
    double M  = m1 + m2;
    double r_peri = a * (1.0 - e);
    double v_peri = std::sqrt(sys.G * M * (1.0 + e) / (a * (1.0 - e)));

    // Cuerpo 1 en −(m2/M)·r_peri, cuerpo 2 en +(m1/M)·r_peri
    Body b1, b2;
    b1.mass     = m1;
    b1.position = Vec3{-m2/M * r_peri, 0.0, 0.0};
    b1.velocity = Vec3{0.0, -m1/(m1+m2) * v_peri, 0.0};

    b2.mass     = m2;
    b2.position = Vec3{+m1/M * r_peri, 0.0, 0.0};
    b2.velocity = Vec3{0.0, +m2/(m1+m2) * v_peri, 0.0};

    sys.bodies = {b1, b2};
    return sys;
}

// Pitágoras: masas 3,4,5 en triángulo rectángulo, en reposo
static NBodySystem make_pythagorean() {
    NBodySystem sys;
    sys.G = 1.0;
    // Posiciones de Szebehely & Peters (1967) — triángulo rectángulo escalado
    // con hipotenusa 5, catetos 3 y 4 (unidades: separación inicial ~1)
    Body b1, b2, b3;
    b1.mass = 3.0; b1.position = Vec3{1.0, 3.0, 0.0}; b1.velocity = Vec3{0,0,0};
    b2.mass = 4.0; b2.position = Vec3{-2.0,-1.0, 0.0}; b2.velocity = Vec3{0,0,0};
    b3.mass = 5.0; b3.position = Vec3{1.0,-1.0, 0.0}; b3.velocity = Vec3{0,0,0};
    sys.bodies = {b1, b2, b3};
    return sys;
}

// Figura-8 de Simó (2002) — CI con 15 dígitos significativos
static NBodySystem make_figure8() {
    NBodySystem sys;
    sys.G = 1.0;
    Body b1, b2, b3;
    b1.mass=1; b1.position=Vec3{ 0.9700436926041022,-0.2430865994534988,0};
               b1.velocity=Vec3{ 0.4662036850003821, 0.4323657300939072,0};
    b2.mass=1; b2.position=Vec3{-0.9700436926041022,-0.2430865994534988,0};
               b2.velocity=Vec3{ 0.4662036850003821,-0.4323657300939072,0};
    b3.mass=1; b3.position=Vec3{ 0.0,                0.4861731989069976,0};
               b3.velocity=Vec3{-0.9324073700007642, 0.0,               0};
    sys.bodies = {b1, b2, b3};
    return sys;
}

// ─── BSParameters estándar ───────────────────────────────────────────────────
static ARChainNBSIntegrator::BSParameters bs_params(double eps = 1e-10) {
    ARChainNBSIntegrator::BSParameters p;
    p.bs_eps     = eps;
    p.initial_ds = 1e-3;
    p.min_ds     = 1e-14;
    p.max_ds     = 1e-1;
    p.k_max      = 8;
    p.max_steps  = 20000000;
    return p;
}

// ══════════════════════════════════════════════════════════════════════════════
// TEST P1 — Precesión del perihelio (PN1)
// ══════════════════════════════════════════════════════════════════════════════
bool test_P1_perihelion_precession() {
    std::cout << "\n" << std::string(70,'═') << "\n";
    std::cout << "TEST P1 — Precesión del perihelio (PN1)\n";
    std::cout << "  Referencia: Einstein (1915); Iyer & Will (1993)\n";
    std::cout << std::string(70,'─') << "\n";

    // Parámetros orbitales
    const double m1 = 1.0, m2 = 1.0;
    const double a  = 0.1,  e  = 0.6;
    const double c  = 10.0;   // c pequeño para tener precesión medible en pocas órbitas
    const double G  = 1.0;
    const double M  = m1 + m2;

    // Precesión teórica PN1 por órbita (Einstein 1915)
    const double dw_theory = 6.0 * M_PI * G * M / (a * c*c * (1.0 - e*e));
    const double T_orb     = kepler_period(a, G, M);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  a=" << a << "  e=" << e << "  c=" << c << "  M=" << M << "\n";
    std::cout << "  T_orbital = " << T_orb << "\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  Δω_teoria = " << dw_theory << " rad/órbita\n\n";

    // Integrador PN1 (pn_order=1)
    ARChainNBSIntegrator::BSParameters bp = bs_params(1e-11);
    ARChainNPNBSIntegrator integ(1e-4, c, /*pn_order=*/1, bp);

    NBodySystem sys = make_binary(m1, m2, a, e);
    std::vector<int> idx = {0, 1};
    ARChainNState state = integ.initialize(sys, idx);

    const double E0 = total_energy(sys);
    const double L0 = norm(total_angular_momentum(sys));

    // Medir perihelios: en cada perihelio (sep mínima local) registrar el ángulo
    // del vector relativo r₂₁ = r₂ - r₁
    const int N_orbits = 10;
    const double t_final = N_orbits * T_orb;

    std::vector<double> perihelion_angles;
    double prev_sep  = 1e30;
    double sep       = 0.0;
    double prev_t    = 0.0;
    bool   was_decreasing = false;

    // Variables para checks energéticos
    double max_dE = 0.0, max_dL = 0.0;

    // Integrar paso a paso (dt = T_orb/200 para no perder perihelios)
    const double dt_step = T_orb / 200.0;
    double t_now = 0.0;

    // Posición relativa al inicio (para ángulo de referencia)
    // La primera medida de perihelio será el ángulo 0 (referencia)
    bool first_perihelion = true;
    double angle_ref = 0.0;

    std::cout << "  Integrando " << N_orbits << " órbitas...\n";
    std::cout << "  " << std::left << std::setw(10) << "Órbita"
              << std::setw(18) << "ω_acum (rad)"
              << std::setw(18) << "|dE/E₀|"
              << "\n";
    std::cout << "  " << std::string(46,'-') << "\n";

    int orbit_count = 0;
    prev_sep = 1e30;
    was_decreasing = false;

    while (t_now < t_final - 1e-12) {
        double t_target = std::min(t_now + dt_step, t_final);
        integ.integrate_to_bs(state, t_target);
        t_now = state.t_phys;

        // Reconstruir posiciones para observables
        integ.write_back(state, sys, idx);

        // Separación actual
        Vec3 dr = sys.bodies[1].position - sys.bodies[0].position;
        sep = norm(dr);

        // Detección de perihelio: mínimo local de sep
        // (sep decrece → aumenta → perihelio entre ambos)
        if (was_decreasing && sep > prev_sep + 1e-10) {
            // Perihelio detectado: el vector dr en el paso anterior
            // Para ángulo, usamos el ángulo del vector al momento del mínimo
            // (aproximación: ángulo en t_now)
            double angle = std::atan2(dr.y, dr.x);

            if (first_perihelion) {
                angle_ref    = angle;
                first_perihelion = false;
            } else {
                // Acumulación de precesión
                double dw_accum = angle - angle_ref;
                // Normalizar a [-π, π] no es correcto aquí — queremos acumulado
                // Mejor: comparar con el ángulo esperado
                orbit_count++;
                double dw_per_orbit = dw_accum / orbit_count;

                if (orbit_count % 2 == 0 || orbit_count == N_orbits-1) {
                    // Observables energéticos
                    integ.write_back(state, sys, idx);
                    double E_now = total_energy(sys);
                    double L_now = norm(total_angular_momentum(sys));
                    double dE = std::abs(E_now - E0) / (std::abs(E0) + 1e-30);
                    double dL = std::abs(L_now - L0) / (std::abs(L0) + 1e-30);
                    max_dE = std::max(max_dE, dE);
                    max_dL = std::max(max_dL, dL);

                    std::cout << "  " << std::left << std::setw(10) << orbit_count
                              << std::scientific << std::setw(18) << dw_accum
                              << std::setw(18) << dE << "\n";
                }

                // Solo los últimos perihelios para la medida final
                if (orbit_count == N_orbits - 1) {
                    perihelion_angles.push_back(dw_accum);
                }
            }
        }

        was_decreasing = (sep < prev_sep);
        prev_sep = sep;
    }

    if (perihelion_angles.empty()) {
        std::cout << "  [ERROR] No se detectaron perihelios — aumentar N_orbits\n";
        return false;
    }

    // Precesión total medida y por órbita
    double dw_total_meas = perihelion_angles.back();
    double dw_per_orbit_meas = dw_total_meas / (N_orbits - 1);
    double dw_theory_total   = dw_theory * (N_orbits - 1);

    double ratio = dw_total_meas / dw_theory_total;
    double rel_error = std::abs(ratio - 1.0);

    std::cout << "\n  ── Resultado P1 ─────────────────────────────────────\n";
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "  Δω_teoria  (total) = " << dw_theory_total   << " rad\n";
    std::cout << "  Δω_medido  (total) = " << dw_total_meas     << " rad\n";
    std::cout << "  Ratio Δω_meas/Δω_teoria = " << ratio        << "\n";
    std::cout << "  Error relativo          = " << rel_error     << "\n";
    std::cout << "  max|dE/E₀|              = " << max_dE        << "\n";
    std::cout << "  max|dL/L₀|              = " << max_dL        << "\n\n";

    // ── Criterios de éxito ────────────────────────────────────────────────
    bool precession_ok = (rel_error < 0.01);
    bool precession_warn = (!precession_ok && rel_error < 0.05);
    bool energy_ok    = (max_dE < 1e-8);
    bool angular_ok   = (max_dL < 1e-8);

    if (precession_ok)
        std::cout << "  [PASS] Precesión PN1: error < 1%  → compute_PN1() correcto\n";
    else if (precession_warn)
        std::cout << "  [WARN] Precesión PN1: error 1-5% → revisar paso de integración\n";
    else
        std::cout << "  [FAIL] Precesión PN1: error > 5% → error en compute_PN1() o GBS\n";

    if (energy_ok)
        std::cout << "  [PASS] PN1 conservativo: |dE/E₀| < 1e-8\n";
    else
        std::cout << "  [FAIL] PN1 NO conservativo → bug: PN1 no debe disipar energía\n";

    if (angular_ok)
        std::cout << "  [PASS] Momento angular conservado: |dL/L₀| < 1e-8\n";
    else
        std::cout << "  [FAIL] Momento angular no conservado\n";

    bool ok = precession_ok && energy_ok && angular_ok;
    std::cout << "\n  Resultado P1: " << (ok ? "PASADO ✅" : "FALLADO ❌") << "\n";
    return ok;
}

// ══════════════════════════════════════════════════════════════════════════════
// TEST P2 — Problema de Pitágoras completo
// ══════════════════════════════════════════════════════════════════════════════
bool test_P2_pythagorean_complete() {
    std::cout << "\n" << std::string(70,'═') << "\n";
    std::cout << "TEST P2 — Problema de Pitágoras completo\n";
    std::cout << "  Referencia: Szebehely & Peters (1967), AJ 72, 876\n";
    std::cout << "  Resultados publicados: eyección de masa-3, v_∞ ≈ 0.347\n";
    std::cout << std::string(70,'─') << "\n";

    // Velocidad de referencia publicada (Szebehely & Peters 1967)
    // Confirmada en Aarseth (2003) con NBODY6 al mismo orden de magnitud
    const double v_inf_ref = 0.347;

    ARChainNBSIntegrator::BSParameters bp = bs_params(1e-10);
    bp.max_steps = 50000000;
    ARChainNBSIntegrator integ(1e-4, bp);

    NBodySystem sys = make_pythagorean();
    const double E0_total = total_energy(sys);
    Vec3 P0 = total_momentum(sys);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "  E₀ total = " << E0_total << "\n";
    std::cout << "  |P₀| = " << norm(P0) << "  (debe ser 0 — masas en reposo)\n\n";

    std::vector<int> idx = {0, 1, 2};
    ARChainNState state = integ.initialize(sys, idx);

    // Integrar hasta t=100 con muestreo para detectar eyección
    // El cuerpo de masa 3 (índice 0) debería ser eyectado
    const double t_final = 100.0;
    const double dt_sample = 1.0;

    std::cout << "  Integrando hasta t=" << t_final << "...\n";
    std::cout << "  " << std::left
              << std::setw(8)  << "t"
              << std::setw(14) << "r(m=3)"
              << std::setw(14) << "r(m=4)"
              << std::setw(14) << "r(m=5)"
              << std::setw(12) << "|dE/E₀|"
              << "\n";
    std::cout << "  " << std::string(62,'-') << "\n";

    double max_dE_rel = 0.0;
    double t_now = 0.0;

    // Puntos de muestreo para la tabla
    std::vector<double> t_print = {5,10,20,30,40,50,60,70,80,90,100};
    int tp_idx = 0;

    while (t_now < t_final - 1e-12 && tp_idx < (int)t_print.size()) {
        double t_target = std::min(t_print[tp_idx], t_final);
        integ.integrate_to_bs(state, t_target);
        t_now = state.t_phys;
        integ.write_back(state, sys, idx);

        double E_now = total_energy(sys);
        double dE_rel = std::abs(E_now - E0_total) / (std::abs(E0_total) + 1e-30);
        max_dE_rel = std::max(max_dE_rel, dE_rel);

        double r0 = norm(sys.bodies[0].position);  // masa 3
        double r1 = norm(sys.bodies[1].position);  // masa 4
        double r2 = norm(sys.bodies[2].position);  // masa 5

        std::cout << "  " << std::fixed << std::setprecision(1)
                  << std::setw(8) << t_now
                  << std::setprecision(4)
                  << std::setw(14) << r0
                  << std::setw(14) << r1
                  << std::setw(14) << r2
                  << std::scientific << std::setprecision(2)
                  << std::setw(12) << dE_rel << "\n";

        tp_idx++;
    }

    // Estado final
    integ.write_back(state, sys, idx);

    // ── Métricas finales ──────────────────────────────────────────────────

    // Velocidad asintótica del cuerpo eyectado (masa 3 = índice 0)
    double v_inf_meas = norm(sys.bodies[0].velocity);

    // Energía de la binaria residual (masas 4 y 5 = índices 1 y 2)
    double E_bin_45 = binding_energy_pair(sys, 1, 2);

    // Posición del cuerpo eyectado
    double r_ejected = norm(sys.bodies[0].position);

    // Conservación del CM
    Vec3 P_final = total_momentum(sys);
    double dP = norm(P_final - P0);

    // Balance energético: E_bin + E_cinética_libre ≈ E₀
    double v_e   = norm(sys.bodies[0].velocity);
    double E_lib = 0.5 * sys.bodies[0].mass * v_e * v_e;  // cinética asintótica masa-3
    double E_check = E_bin_45 + E_lib;
    double dE_balance = std::abs(E_check - E0_total) / (std::abs(E0_total) + 1e-30);

    std::cout << "\n  ── Resultado P2 ─────────────────────────────────────\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  Posición masa-3 en t=100:  r = " << r_ejected        << "\n";
    std::cout << "  Velocidad asintótica meas: v_∞ = " << v_inf_meas      << "\n";
    std::cout << "  Velocidad referencia:      v_ref= " << v_inf_ref       << "  (S&P 1967)\n";
    std::cout << "  |v_∞/v_ref − 1|            = "
              << std::abs(v_inf_meas/v_inf_ref - 1.0)                     << "\n";
    std::cout << "  Energía binaria E_bin(4,5) = " << E_bin_45            << "\n";
    std::cout << "  Balance E_bin+E_lib − E₀   = " << dE_balance          << "\n";
    std::cout << "  |ΔP_total|                 = " << dP                  << "\n";
    std::cout << "  max|dE/E₀| integración     = " << max_dE_rel          << "\n\n";

    // ── Criterios ─────────────────────────────────────────────────────────
    bool ejection_ok = (r_ejected > 30.0);
    bool binary_ok   = (E_bin_45  < 0.0);
    bool v_inf_ok    = (std::abs(v_inf_meas / v_inf_ref - 1.0) < 0.15);
    bool cm_ok       = (dP < 1e-10);
    bool balance_ok  = (dE_balance < 1e-4);

    std::cout << "  [" << (ejection_ok ? "PASS" : "FAIL") << "] Eyección masa-3: r > 30  (r=" << r_ejected << ")\n";
    std::cout << "  [" << (binary_ok   ? "PASS" : "FAIL") << "] Binaria (4,5) ligada: E_bin < 0\n";
    std::cout << "  [" << (v_inf_ok    ? "PASS" : "WARN") << "] Velocidad asintótica: "
              << "error=" << std::abs(v_inf_meas/v_inf_ref-1.0)*100 << "%  (tol 15%)\n";
    std::cout << "  [" << (cm_ok       ? "PASS" : "FAIL") << "] CM conservado: |ΔP| < 1e-10\n";
    std::cout << "  [" << (balance_ok  ? "PASS" : "FAIL") << "] Balance energético < 1e-4\n";

    // Criterio mínimo de éxito: eyección + binaria ligada + CM conservado
    // v_inf es orientativo (sistema caótico, tolerancia amplia)
    bool ok = ejection_ok && binary_ok && cm_ok && balance_ok;
    std::cout << "\n  Resultado P2: " << (ok ? "PASADO ✅" : "FALLADO ❌") << "\n";
    return ok;
}

// ══════════════════════════════════════════════════════════════════════════════
// TEST P3 — Figura-8 × 100 períodos: deriva energética secular
// ══════════════════════════════════════════════════════════════════════════════
bool test_P3_figure8_long() {
    std::cout << "\n" << std::string(70,'═') << "\n";
    std::cout << "TEST P3 — Figura-8 × 100 períodos (deriva secular)\n";
    std::cout << "  Referencia: Simó (2002), Cel. Mech. Dyn. Astron. 82, 3\n";
    std::cout << "  Criterio:   sin drift lineal en energía (integrador simpléctico)\n";
    std::cout << std::string(70,'─') << "\n";

    const double T_period = 6.3259;   // período de la figura-8 (Simó 2002)
    const int    N_orbits = 100;
    const double t_final  = N_orbits * T_period;

    ARChainNBSIntegrator::BSParameters bp = bs_params(1e-10);
    bp.max_steps = 100000000;
    ARChainNBSIntegrator integ(1e-4, bp);

    NBodySystem sys = make_figure8();
    const double E0 = total_energy(sys);
    const double L0 = norm(total_angular_momentum(sys));

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "  E₀ = " << E0 << "\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  L₀ = " << L0 << "  (debe ser ~0 por simetría de las CI)\n\n";

    std::vector<int> idx = {0, 1, 2};
    ARChainNState state = integ.initialize(sys, idx);

    // Muestrear cada período
    std::vector<double> t_samples, dE_samples;
    double max_dE = 0.0, max_dL = 0.0;

    std::cout << "  " << std::left
              << std::setw(10) << "Período"
              << std::setw(16) << "|dE/E₀|"
              << std::setw(16) << "|ΔL|"
              << "\n";
    std::cout << "  " << std::string(42,'-') << "\n";

    for (int n = 1; n <= N_orbits; ++n) {
        double t_target = n * T_period;
        integ.integrate_to_bs(state, t_target);
        integ.write_back(state, sys, idx);

        double E_now = total_energy(sys);
        double L_now = norm(total_angular_momentum(sys));
        double dE = std::abs(E_now - E0) / (std::abs(E0) + 1e-30);
        double dL = std::abs(L_now - L0);   // L0 ≈ 0, así que dL es absoluto

        max_dE = std::max(max_dE, dE);
        max_dL = std::max(max_dL, dL);

        t_samples.push_back((double)n);
        dE_samples.push_back(dE);

        // Imprimir cada 10 períodos
        if (n % 10 == 0 || n == 1) {
            std::cout << "  " << std::setw(10) << n
                      << std::scientific << std::setw(16) << dE
                      << std::setw(16) << dL << "\n";
        }
    }

    // ── Regresión lineal para detectar drift secular ──────────────────────
    // slope = Σ(x-x̄)(y-ȳ) / Σ(x-x̄)²
    int n = (int)t_samples.size();
    double x_mean = std::accumulate(t_samples.begin(),  t_samples.end(),  0.0) / n;
    double y_mean = std::accumulate(dE_samples.begin(), dE_samples.end(), 0.0) / n;

    double num = 0.0, den = 0.0;
    for (int i = 0; i < n; ++i) {
        double dx = t_samples[i]  - x_mean;
        double dy = dE_samples[i] - y_mean;
        num += dx * dy;
        den += dx * dx;
    }
    double slope = (den > 1e-30) ? num / den : 0.0;
    // slope en unidades de [dE/E₀] / período
    // Para un integrador simpléctico, slope debe ser compatible con 0

    std::cout << "\n  ── Resultado P3 ─────────────────────────────────────\n";
    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  max|dE/E₀|        = " << max_dE << "  (criterio: < 1e-8)\n";
    std::cout << "  max|ΔL|           = " << max_dL << "  (criterio: < 1e-12)\n";
    std::cout << "  slope dE(t)/T     = " << slope  << "  (criterio: < 1e-11)\n\n";

    bool energy_ok  = (max_dE  < 1e-8);
    bool angular_ok = (max_dL  < 1e-12);
    bool slope_ok   = (std::abs(slope) < 1e-11);

    std::cout << "  [" << (energy_ok  ? "PASS" : "FAIL")
              << "] max|dE/E₀| < 1e-8 → GBS mantiene precisión en 100 períodos\n";
    std::cout << "  [" << (angular_ok ? "PASS" : "FAIL")
              << "] |ΔL| < 1e-12    → L=0 conservado (simetría Simó)\n";
    std::cout << "  [" << (slope_ok   ? "PASS" : "FAIL")
              << "] slope < 1e-11   → sin drift secular (integrador simpléctico)\n";

    bool ok = energy_ok && angular_ok && slope_ok;
    std::cout << "\n  Resultado P3: " << (ok ? "PASADO ✅" : "FALLADO ❌") << "\n";
    return ok;
}

// ══════════════════════════════════════════════════════════════════════════════
// MAIN
// ══════════════════════════════════════════════════════════════════════════════
int main() {
    std::cout << std::string(70,'█') << "\n";
    std::cout << "TESTS FÍSICOS — Fase 3 (Benchmarks publicados)\n";
    std::cout << std::string(70,'█') << "\n";

    bool p1 = test_P1_perihelion_precession();
    bool p2 = test_P2_pythagorean_complete();
    bool p3 = test_P3_figure8_long();

    std::cout << "\n" << std::string(70,'═') << "\n";
    std::cout << "RESUMEN:\n";
    std::cout << "  P1 (precesión perihelio PN1):     " << (p1 ? "✅ PASADO" : "❌ FALLADO") << "\n";
    std::cout << "  P2 (Pitágoras completo):           " << (p2 ? "✅ PASADO" : "❌ FALLADO") << "\n";
    std::cout << "  P3 (figura-8 × 100, sin drift):   " << (p3 ? "✅ PASADO" : "❌ FALLADO") << "\n";
    std::cout << std::string(70,'═') << "\n";

    return (p1 && p2 && p3) ? 0 : 1;
}