// Diagnostico: violacion del constraint bilineal Q.P y T_mikkola
#include <iostream>
#include <iomanip>
#include <cmath>
#include "nbody_system.h"
#include "chain3_bs_integrator.h"
#include "initial_conditions.h"
#include "chain_ks.h"

int main() {
    std::cout << std::setprecision(15);
    NBodySystem sys0 = InitialConditions::figure_eight();
    double E0 = sys0.total_energy();

    Chain3BSIntegrator::BSParameters params;
    params.eps=1e-10; params.max_steps=1000000; params.verbose=false;
    Chain3BSIntegrator bs(params);
    Chain3State state = bs.initialize(sys0, 0, 1, 2);

    // Mediciones en t=0
    auto measure = [&](const char* label, const Chain3State& s) {
        double m1=s.m1(), m2=s.m2(), m3=s.m3();
        Vec3 W1=ks_to_w(s.Q1,s.P1), W2=ks_to_w(s.Q2,s.P2);
        double T11=1./m1+1./m2, T12=-1./m2, T22=1./m2+1./m3;
        Vec3 A1=T11*W1+T12*W2, A2=T12*W1+T22*W2;
        double T_mik = A1.dot(W1) + A2.dot(W2);
        Vec3 R1=ks_to_r(s.Q1), R2=ks_to_r(s.Q2), R3=R1+R2;
        double U = -(m1*m2/R1.norm()+m2*m3/R2.norm()+m1*m3/R3.norm());
        // Bilinear constraint violation
        double qp1 = s.Q1.x*s.P1.x+s.Q1.y*s.P1.y+s.Q1.z*s.P1.z+s.Q1.w*s.P1.w;
        double qp2 = s.Q2.x*s.P2.x+s.Q2.y*s.P2.y+s.Q2.z*s.P2.z+s.Q2.w*s.P2.w;
        double q1n = s.Q1.norm2(), q2n = s.Q2.norm2();
        std::cout << label << ":\n";
        std::cout << "  T_mik=" << T_mik << "  U=" << U
                  << "  H_mik=" << T_mik+U << "  state.energy=" << s.energy << "\n";
        std::cout << "  E_phys=0.5*T_mik+U=" << 0.5*T_mik+U << "  (E0=" << E0 << ")\n";
        std::cout << "  Q1.P1/|Q1|^2=" << qp1/q1n << "  Q2.P2/|Q2|^2=" << qp2/q2n
                  << "  (should be 0)\n";
        std::cout << "  |Q1|^2=" << q1n << "  |Q2|^2=" << q2n << "\n\n";
    };

    measure("t=0", state);

    // Integrar pasos pequeños y medir
    double t_targets[] = {0.001, 0.01, 0.1, 1.0, 6.325919634};
    for (double T : t_targets) {
        double t_done;
        bs.integrate_precise(state, T, t_done, sys0);
        char label[64];
        snprintf(label, sizeof(label), "t=%.4f", t_done);
        measure(label, state);
    }

    return 0;
}