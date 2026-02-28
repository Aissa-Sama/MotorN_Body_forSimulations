// regularization/chain/chain3_bs_integrator.cpp
#include "chain3_bs_integrator.h"
#include "chain_ks.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

// ============================================================================
// CONSTRUCTOR
// ============================================================================
Chain3BSIntegrator::Chain3BSIntegrator(const BSParameters& params)
    : Chain3Integrator(1e-4)   // fallback dtau base (raramente usado)
    , bs_params_(params)
{}

// NOTA sobre el constraint bilineal Q·P = 0:
// En la cadena de Mikkola, se puede demostrar que Q·P = 2·R·W
// donde R es la separación física del eslabón y W el momento conjugado.
// Este valor NO es cero en general — es un INVARIANTE DE GAUGE, no una condición
// dinámica. Las ecuaciones de movimiento dQ/dτ y dP/dτ son independientes de Q·P.
//
// Proyectar P → P - (Q·P/|Q|²)·Q cambiaría W según:
//   ΔW = -(Q·P / |Q|⁴) · R
// introduciendo un error físico real en la velocidad y en T_mikkola.
// Por esta razón NO se aplica ninguna proyección del constraint bilineal.

// ============================================================================
// PASO MEDIO MODIFICADO (Gragg / Bulirsch-Stoer)
//
// Implementa el "modified midpoint method":
//   z0 = y0
//   z1 = z0 + h * f(z0)                    ← primer paso Euler
//   z_{m+1} = z_{m-1} + 2h * f(z_m)        ← pasos medios
//   result = (z_n + z_{n-1} + h * f(z_n)) / 2  ← paso final suavizante
//
// IMPORTANTE: cm_time y tau se propagan exactamente igual que en rk4_step,
// a través de deriv.cm_time = g (el factor de regularización).
// NO se asigna manualmente cm_time = y0.cm_time + dt.
// ============================================================================
void Chain3BSIntegrator::modified_midpoint_step(
    const Chain3State& y0,
    double dtau,
    int n_steps,
    const NBodySystem& full_system,
    Chain3State& result
) {
    Chain3State z_prev = y0;
    Chain3State z_curr = y0;

    // --- Primer paso (Euler) ---
    Chain3State d0;
    compute_derivatives(z_prev, full_system, d0);

    z_curr.Q1      = z_prev.Q1      + dtau * d0.Q1;
    z_curr.P1      = z_prev.P1      + dtau * d0.P1;
    z_curr.Q2      = z_prev.Q2      + dtau * d0.Q2;
    z_curr.P2      = z_prev.P2      + dtau * d0.P2;
    // cm_time avanza con g*dtau (derivada = g)
    double dcm0    = dtau * d0.cm_time;
    z_curr.cm_time = z_prev.cm_time + dcm0;
    // cm_pos avanza con cm_vel * dt_fisico (igual que RK4)
    z_curr.cm_pos  = z_prev.cm_pos  + z_prev.cm_vel * dcm0;
    z_curr.tau     = z_prev.tau     + dtau;

    // --- Pasos medios ---
    for (int m = 1; m < n_steps; ++m) {
        Chain3State dc;
        compute_derivatives(z_curr, full_system, dc);

        Chain3State z_next = z_curr;
        z_next.Q1      = z_prev.Q1      + 2.0 * dtau * dc.Q1;
        z_next.P1      = z_prev.P1      + 2.0 * dtau * dc.P1;
        z_next.Q2      = z_prev.Q2      + 2.0 * dtau * dc.Q2;
        z_next.P2      = z_prev.P2      + 2.0 * dtau * dc.P2;
        double dcm     = 2.0 * dtau * dc.cm_time;
        z_next.cm_time = z_prev.cm_time + dcm;
        z_next.cm_pos  = z_prev.cm_pos  + z_curr.cm_vel * dcm;
        z_next.tau     = z_curr.tau     + dtau;

        z_prev = z_curr;
        z_curr = z_next;
    }

    // --- Paso final suavizante ---
    Chain3State dn;
    compute_derivatives(z_curr, full_system, dn);

    double dcmn    = dtau * dn.cm_time;
    result.Q1      = 0.5 * (z_curr.Q1      + z_prev.Q1      + dtau * dn.Q1);
    result.P1      = 0.5 * (z_curr.P1      + z_prev.P1      + dtau * dn.P1);
    result.Q2      = 0.5 * (z_curr.Q2      + z_prev.Q2      + dtau * dn.Q2);
    result.P2      = 0.5 * (z_curr.P2      + z_prev.P2      + dtau * dn.P2);
    result.cm_time = 0.5 * (z_curr.cm_time + z_prev.cm_time + dcmn);
    result.cm_pos  = 0.5 * (z_curr.cm_pos  + z_prev.cm_pos  + y0.cm_vel * dcmn);  // CORRECCIÓN: usar y0.cm_vel (cm_vel no varía durante el midpoint, pero y0.cm_vel es más claro)
    result.cm_vel  = y0.cm_vel;
    result.tau     = y0.tau + (double)n_steps * dtau;
    result.energy    = y0.energy;
    result.masses[0] = y0.masses[0];
    result.masses[1] = y0.masses[1];
    result.masses[2] = y0.masses[2];

}

// ============================================================================
// EXTRAPOLACIÓN POLINÓMICA DE RICHARDSON (Neville-Aitken)
//
// T[k][m] = T[k][m-1] + (T[k][m-1] - T[k-1][m-1]) / ((h_{k-m}/h_k)^2 - 1)
// ============================================================================
void Chain3BSIntegrator::polynomial_extrapolation(
    const std::vector<Chain3State>& results,
    const std::vector<double>&      step_sizes,
    int                             k_current,
    Chain3State&                    extrapolated
) const {
    // Copia de trabajo: T[j] = resultado del nivel j
    std::vector<Chain3State> T(results.begin(),
                                results.begin() + k_current + 1);

    for (int m = 1; m <= k_current; ++m) {
        for (int k = k_current; k >= m; --k) {
            double ratio  = step_sizes[k - m] / step_sizes[k];
            double factor = ratio * ratio - 1.0;
            if (std::abs(factor) < 1e-15) factor = 1e-15;

            T[k].Q1      = T[k].Q1      + (T[k].Q1      - T[k-1].Q1)      / factor;
            T[k].P1      = T[k].P1      + (T[k].P1      - T[k-1].P1)      / factor;
            T[k].Q2      = T[k].Q2      + (T[k].Q2      - T[k-1].Q2)      / factor;
            T[k].P2      = T[k].P2      + (T[k].P2      - T[k-1].P2)      / factor;
            T[k].cm_time = T[k].cm_time + (T[k].cm_time - T[k-1].cm_time) / factor;
        }
    }

    extrapolated = T[k_current];
}

// ============================================================================
// ESTIMACIÓN DE ERROR ENTRE DOS NIVELES DE EXTRAPOLACIÓN
// ============================================================================
double Chain3BSIntegrator::estimate_error_bs(
    const Chain3State& high_order,
    const Chain3State& low_order
) const {
    double scale_Q1 = std::max(high_order.Q1.norm(), 1e-10);
    double scale_P1 = std::max(high_order.P1.norm(), 1e-10);
    double scale_Q2 = std::max(high_order.Q2.norm(), 1e-10);
    double scale_P2 = std::max(high_order.P2.norm(), 1e-10);

    return ((high_order.Q1 - low_order.Q1).norm() / scale_Q1
          + (high_order.P1 - low_order.P1).norm() / scale_P1
          + (high_order.Q2 - low_order.Q2).norm() / scale_Q2
          + (high_order.P2 - low_order.P2).norm() / scale_P2) / 4.0;
}

// ============================================================================
// FACTOR DE ESCALA PARA PASO NUEVO
// ============================================================================
double Chain3BSIntegrator::step_scale_factor(double error, int k_opt) const {
    if (error < 1e-30) return 4.0;  // error ~0 → máximo crecimiento
    double exponent = 1.0 / (2.0 * k_opt + 1.0);
    double scale    = bs_params_.safety * std::pow(1.0 / error, exponent);
    return std::clamp(scale, 0.1, 4.0);
}

// ============================================================================
// CÁLCULO DE ENERGÍA DESDE ESTADO (para diagnóstico)
// ============================================================================
double Chain3BSIntegrator::compute_energy_from_state(const Chain3State& s) const {
    Vec3 R1 = ks_to_r(s.Q1), R2 = ks_to_r(s.Q2), R3 = R1 + R2;
    Vec3 W1 = ks_to_w(s.Q1, s.P1), W2 = ks_to_w(s.Q2, s.P2);
    double m1 = s.m1(), m2 = s.m2(), m3 = s.m3();
    double T11 = 1.0/m1 + 1.0/m2, T12 = -1.0/m2, T22 = 1.0/m2 + 1.0/m3;
    Vec3 A1 = T11*W1 + T12*W2, A2 = T12*W1 + T22*W2;
    double T_kin = dot(A1, W1) + dot(A2, W2);
    double U = -(m1*m2/R1.norm() + m2*m3/R2.norm() + m1*m3/R3.norm());
    return T_kin + U;
}


// CONVENCIÓN: Esta función calcula la energía con la cinética clásica T_clásica = ½·T_Mikkola,
// que equivale a Σ½·mᵢ·vᵢ² en coordenadas físicas. Es DIFERENTE a state.energy (convención
// Mikkola) que usa T_Mikkola sin factor ½. No usar para comparar con state.energy.
double Chain3BSIntegrator::compute_classical_energy(const Chain3State& s) const {
    Vec3 R1 = ks_to_r(s.Q1), R2 = ks_to_r(s.Q2), R3 = R1 + R2;
    Vec3 W1 = ks_to_w(s.Q1, s.P1), W2 = ks_to_w(s.Q2, s.P2);
    double m1 = s.m1(), m2 = s.m2(), m3 = s.m3();
    double T11 = 1.0/m1 + 1.0/m2, T12 = -1.0/m2, T22 = 1.0/m2 + 1.0/m3;
    Vec3 A1 = T11*W1 + T12*W2, A2 = T12*W1 + T22*W2;
    double T_mikkola = dot(A1, W1) + dot(A2, W2);
    double U = -(m1*m2/R1.norm() + m2*m3/R2.norm() + m1*m3/R3.norm());
    return 0.5 * T_mikkola + U;  // energía clásica: T_clásica = ½·T_Mikkola
}

double Chain3BSIntegrator::compute_energy_error(
    const Chain3State& state, double E_ref) const
{
    return std::abs(compute_energy_from_state(state) - E_ref) / std::abs(E_ref);
}

// ============================================================================
// UN PASO BS COMPLETO
// Error = extrap[k] - extrap[k-1]. No se proyecta Q·P (ver nota arriba).
// ============================================================================
bool Chain3BSIntegrator::bs_step(
    Chain3State& state,
    double       dtau_total,
    double&      error_estimate,
    double&      dtau_out,
    const NBodySystem& full_system
) {
    const Chain3State y_save = state;

    std::vector<Chain3State> raw;
    std::vector<double>      step_sizes;
    raw.reserve(K_MAX);
    step_sizes.reserve(K_MAX);

    Chain3State extrap_prev;
    bool has_prev  = false;
    error_estimate = 1e30;

    for (int k = 0; k < K_MAX; ++k) {
        int    n      = n_seq_[k];
        double dtau_k = dtau_total / n;

        Chain3State r;
        try {
            modified_midpoint_step(y_save, dtau_k, n, full_system, r);
        } catch (const std::exception&) {
            dtau_out = dtau_total * 0.25;
            return false;
        }
        raw.push_back(r);
        step_sizes.push_back(dtau_k);

        if (k >= 1) {
            Chain3State extrap_curr;
            polynomial_extrapolation(raw, step_sizes, k, extrap_curr);

            if (has_prev) {
                error_estimate = estimate_error_bs(extrap_curr, extrap_prev);
                double scaled  = error_estimate / bs_params_.eps;

                if (scaled <= 1.0) {
                    state = extrap_curr;
                    dtau_out = dtau_total * step_scale_factor(scaled, k);
                    return true;
                }
                if (scaled > 1000.0 && k >= 3) break;
            }
            extrap_prev = extrap_curr;
            has_prev    = true;
        }
    }

    // No convergió: aceptar con paso reducido
    if (has_prev) {
        state = extrap_prev;
        dtau_out = dtau_total * 0.5;
        return true;
    }
    dtau_out = dtau_total * 0.5;
    return false;
}

// ============================================================================
// INTEGRACIÓN PRINCIPAL
//
// El control de paso vive aquí, en tiempo ficticio (dtau).
// La relación dt_phys ≈ g * dtau donde g = 1/L se integra automáticamente
// a través de deriv.cm_time en modified_midpoint_step.
// ============================================================================
void Chain3BSIntegrator::integrate_precise(
    Chain3State&       state,
    double             t_target,
    double&            t_achieved,
    const NBodySystem& full_system
) {
    const double E0_mikkola  = state.energy;
    double dtau    = bs_params_.initial_dtau;
    double max_dE_ = 0.0;

    int total_steps              = 0;
    int rejected_steps           = 0;
    int consecutive_singularities = 0;

    while (state.cm_time < t_target - 1e-14) {
        double t_left = t_target - state.cm_time;
        if (t_left < 1e-14) break;

        double error   = 0.0;
        double dtau_out = dtau;
        Chain3State trial = state;

        bool accepted = false;
        try {
            accepted = bs_step(trial, dtau, error, dtau_out, full_system);
            consecutive_singularities = 0;
        } catch (const std::runtime_error&) {
            // Singularidad numerica: reducir paso y reintentar
            dtau *= 0.25;
            dtau  = std::max(dtau, bs_params_.min_dtau);
            ++consecutive_singularities;
            ++rejected_steps;
            if (consecutive_singularities > 20)
                throw std::runtime_error("Chain3BS: singularidad persistente");
            continue;
        }

        if (accepted) {
            if (trial.cm_time > t_target + 1e-12) {
                dtau = std::max(dtau * 0.5, bs_params_.min_dtau);
                continue;
            }
            // Control de energía: si |ΔH/H₀| > energy_tol, reducir paso y reintentar.
            // Esto previene que la deriva de Q·P (Q·P = 2·R·W ≠ 0) se amplifique
            // exponencialmente en la extrapolación de Richardson.
            // IMPORTANTE: el rechazo solo opera si aún hay margen de paso Y de intentos.
            double E_trial = compute_energy_from_state(trial);
            double dE_trial = std::abs(E_trial - E0_mikkola) / (std::abs(E0_mikkola) + 1e-30);
            if (dE_trial > bs_params_.energy_tol
                && dtau > bs_params_.min_dtau * 4.0
                && rejected_steps < bs_params_.max_steps / 2) {
                dtau = std::max(dtau * 0.5, bs_params_.min_dtau);
                ++rejected_steps;
                continue;
            }
            state = trial;
            dtau  = std::clamp(dtau_out, bs_params_.min_dtau, bs_params_.max_dtau);
            ++total_steps;
            if (dE_trial > max_dE_) max_dE_ = dE_trial;
        } else {
            dtau = std::clamp(dtau_out, bs_params_.min_dtau, bs_params_.max_dtau);
            ++rejected_steps;
        }

        if (total_steps > bs_params_.max_steps)
            throw std::runtime_error("Chain3BS: maximo de pasos excedido");
    }

    t_achieved = state.cm_time;

    if (bs_params_.verbose) {
        std::cout << "[BS] " << total_steps << " pasos, "
                  << rejected_steps << " rechazados"
                  << "  max|dE/E| = " << max_dE_
                  << "\n";
    }
}