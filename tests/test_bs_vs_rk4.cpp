// tests/test_bs_vs_rk4.cpp
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include "nbody_system.h"
#include "chain3_bs_integrator.h"
#include "chain3_integrator.h"

// Binaria circular exacta: m1=m2=1, separacion=1
// v_orbital = sqrt(G*M/4/r) = sqrt(2/4/0.5) = 1.000000
// Periodo = 2*pi*sqrt(a^3/GM) = 4.442883
static NBodySystem make_binary() {
    NBodySystem sys;
    sys.bodies.resize(3);
    sys.bodies[0].mass = 1.0;
    sys.bodies[0].position = {-0.5, 0.0, 0.0};
    sys.bodies[0].velocity = {0.0, -0.7071067812, 0.0};
    sys.bodies[1].mass = 1.0;
    sys.bodies[1].position = { 0.5, 0.0, 0.0};
    sys.bodies[1].velocity = {0.0,  0.7071067812, 0.0};
    // Tercer cuerpo muy lejano y ligero
    sys.bodies[2].mass = 1e-10;
    sys.bodies[2].position = {1000.0, 0.0, 0.0};
    sys.bodies[2].velocity = {0.0, 0.0, 0.0};
    return sys;
}

int main() {
    std::cout << std::setprecision(12);
    const double T_orbit = 4.4428829382;  // 1 periodo
    const double T_total = 10.0 * T_orbit;  // 10 periodos

    std::cout << "=== Conservacion de energia: binaria kepleriana ===\n";
    std::cout << "Periodo orbital: " << T_orbit << "\n";
    std::cout << "Integrando " << T_total << " unidades de tiempo\n\n";

    // ---- RK4 ----
    {
        NBodySystem sys = make_binary();
        double E0 = sys.total_energy();
        std::cout << "E0 = " << E0 << "\n";

        Chain3Integrator rk4(1e-4);
        Chain3State s = rk4.initialize(sys, 0, 1, 2);
        double t_done;
        IntegrationParams p;
        p.abs_tol = 1e-12; p.min_dtau = 1e-10; p.max_dtau = 1e-2;
        rk4.integrate(s, T_total, t_done, p, sys);

        NBodySystem out; out.bodies.resize(3);
        out.bodies[0].mass = s.m1(); out.bodies[1].mass = s.m2(); out.bodies[2].mass = s.m3();
        rk4.write_back(s, out, 0, 1, 2);
        double dE = std::abs(out.total_energy() - E0) / std::abs(E0);
        std::cout << "RK4 |dE/E0| = " << dE << "  t=" << t_done << "\n";
    }

    // ---- BS ----
    {
        NBodySystem sys = make_binary();
        double E0 = sys.total_energy();

        Chain3BSIntegrator::BSParameters p;
        p.eps = 1e-10; p.max_steps = 200000; p.verbose = false;
        Chain3BSIntegrator bs(p);
        Chain3State s = bs.initialize(sys, 0, 1, 2);
        double t_done;
        try {
            bs.integrate_precise(s, T_total, t_done, sys);
        } catch (const std::exception& e) {
            std::cout << "BS exception: " << e.what() << " at t=" << s.cm_time << "\n";
            return 1;
        }

        NBodySystem out; out.bodies.resize(3);
        out.bodies[0].mass = s.m1(); out.bodies[1].mass = s.m2(); out.bodies[2].mass = s.m3();
        bs.write_back(s, out, 0, 1, 2);
        double dE = std::abs(out.total_energy() - E0) / std::abs(E0);
        std::cout << "BS  |dE/E0| = " << dE << "  t=" << t_done << "\n";
    }

    return 0;
}