// regularization/ks/ks_integrator.cpp
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "ks_integrator.h"
#include "ks_transform.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

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
    
    for (int i = 0; i < 4; ++i) {
        ks.u[i] = 0.0;
        ks.w[i] = 0.0;
    }
    
    double m1 = b.mass1();
    double m2 = b.mass2();
    double mu = b.reduced_mass();
    Vec3 r = b.relative_position();
    Vec3 v = b.relative_velocity();
    
    double rnorm = norm(r);
    if (rnorm < 1e-12) {
        throw std::runtime_error("KSIntegrator::to_ks: posición relativa demasiado pequeña");
    }
    
    r_to_ks(r, ks.u);
    v_to_ks(v, ks.u, ks.w);
    
    // ENERGÍA: h = energía específica (completa)
    double h = 0.5 * mu * dot(v, v) - (m1 * m2) / rnorm;   // h = -0.5
    ks.energy = -h;  // energía relativa positiva para órbitas ligadas
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
// UN PASO EN TIEMPO FICTICIO
// ============================================================================
void KSIntegrator::step_ks(KSState& ks) {
    // Ecuaciones: du/dτ = w, dw/dτ = -(energy/2) * u
    
    // Drift
    for (int i = 0; i < 4; ++i) {
        ks.u[i] += internal_dt * ks.w[i];
    }
    
    // Kick - CON FACTOR 1/2
    for (int i = 0; i < 4; ++i) {
        ks.w[i] -= internal_dt * ks.energy * ks.u[i];
    }
    
    ks.tau += internal_dt;
}

// ============================================================================
// INTEGRACIÓN PRINCIPAL
// ============================================================================
void KSIntegrator::integrate(BinaryState& binary, double dt_global) {
    if (dt_global <= 0.0) return;

    KSState ks = to_ks(binary);
    double r = norm(binary.relative_position());
    if (r < 1e-12) {
        throw std::runtime_error("Separación nula en KS");
    }

    // Relación: dτ = dt / (4 r)   →   Δτ_total ≈ dt_global / (4 r)
    // Usamos r inicial como aproximación (válido para órbitas casi circulares)
    double dtau_total = dt_global / r;

    // Número de sub-pasos: queremos que dtau sea pequeño
    // Opción 1: basado en internal_dt (recomendado)
    int n_steps = static_cast<int>(std::ceil(dtau_total / internal_dt));
    if (n_steps < 10) n_steps = 10;  // mínimo razonable

    double dtau = dtau_total / static_cast<double>(n_steps);

    // Bucle de integración (drift-kick simple)
    for (int i = 0; i < n_steps; ++i) {
        // drift (actualiza u)
        for (int j = 0; j < 4; ++j) {
            ks.u[j] += dtau * ks.w[j];
        }
        // kick (actualiza w)
        for (int j = 0; j < 4; ++j) {
            ks.w[j] -= dtau * ks.energy * ks.u[j];
        }
        ks.tau += dtau;
    }
    
    double u_norm2 = 0.0;
    for (int j = 0; j < 4; ++j) u_norm2 += ks.u[j] * ks.u[j];
    double r_final_ks = u_norm2;  // debería ser ≈ r inicial
    std::cout << "  Norma |u|^2 inicial vs final: " << r << " → " << r_final_ks << "\n";
    
    from_ks(ks, binary);
    binary.advance_cm(dt_global);
}