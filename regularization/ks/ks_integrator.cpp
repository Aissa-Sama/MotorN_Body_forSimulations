// regularization/ks/ks_integrator.cpp
#include "ks_integrator.h"
#include "ks_transform.h"
#include "binary_state.h"  // ← ESTE INCLUDE FALTABA
#include <cmath>
#include <algorithm>
#include <stdexcept>

// ============================================================================
// CONSTRUCTOR
// ============================================================================
KSIntegrator::KSIntegrator(double internal_dt)
    : internal_dt(internal_dt) {}

// ============================================================================
// TRANSFORMACIÓN A COORDENADAS KS
// ============================================================================
KSState KSIntegrator::to_ks(const BinaryState& b) const {
    KSState ks;
    
    double m1 = b.mass1();
    double m2 = b.mass2();
    double mu = b.reduced_mass();
    Vec3 r = b.relative_position();
    Vec3 v = b.relative_velocity();

    double rnorm = norm(r);
    if (rnorm < 1e-12)
        throw std::runtime_error("KSIntegrator::to_ks: posición relativa demasiado pequeña");

    r_to_ks(r, ks.u);
    v_to_ks(v, ks.u, ks.w);

    // Energía física (negativa para órbitas ligadas)
    double E = 0.5 * mu * dot(v, v) - (m1 * m2) / rnorm;
    ks.energy = -E;
    ks.tau = 0.0;

    return ks;
}

// ============================================================================
// TRANSFORMACIÓN DESDE COORDENADAS KS
// ============================================================================
void KSIntegrator::from_ks(const KSState& ks, BinaryState& b) const {
    Vec3 new_r = ks_to_r(ks.u);
    Vec3 new_v = ks_to_v(ks.u, ks.w);
    b.set_relative_position(new_r);
    b.set_relative_velocity(new_v);
}

// ============================================================================
// INTEGRACIÓN PRINCIPAL
// ============================================================================
void KSIntegrator::integrate(BinaryState& binary, double dt_global) {
    if (dt_global <= 0.0) return;

    KSState ks = to_ks(binary);

    double r0 = norm(binary.relative_position());
    if (r0 < 1e-12)
        throw std::runtime_error("KSIntegrator::integrate: separación nula");

    // Relación tiempo físico ↔ ficticio: dt = 2|u|² dτ
    double dtau_total = dt_global / (2.0 * r0);

    int n_steps = static_cast<int>(std::ceil(dtau_total / internal_dt));
    if (n_steps < 10) n_steps = 10;

    double dtau = dtau_total / static_cast<double>(n_steps);

    for (int i = 0; i < n_steps; ++i) {
        // Leapfrog en tiempo ficticio
        for (int j = 0; j < 4; ++j)
            ks.w[j] -= 0.5 * dtau * ks.energy * ks.u[j];

        for (int j = 0; j < 4; ++j)
            ks.u[j] += dtau * ks.w[j];

        for (int j = 0; j < 4; ++j)
            ks.w[j] -= 0.5 * dtau * ks.energy * ks.u[j];

        ks.tau += dtau;
    }

    from_ks(ks, binary);
    binary.advance_cm(dt_global);
}