// regularization/chain/chain3_bs_integrator.h
#pragma once
#include "chain3_integrator.h"
#include "chain_ks.h"
#include <vector>
#include <array>
#include <iostream>

// ============================================================================
// INTEGRADOR BULIRSCH-STOER PARA CHAIN3
//
// Reemplaza RK4 por el método de paso medio modificado (Gragg) con
// extrapolación polinómica de Richardson (Neville-Aitken).
//
// Garantiza |ΔE/E| < eps para la integración de la figura-8 y órbitas
// gravitacionales cercanas.
//
// Referencia: Bulirsch & Stoer (1966), Press et al. Numerical Recipes cap. 17
// ============================================================================
class Chain3BSIntegrator : public Chain3Integrator {
public:

    struct BSParameters {
        double eps      = 1e-12;    ///< Tolerancia de error objetivo (extrapolación)
        int    max_steps = 10000;   ///< Máximo de pasos BS totales
        int    max_extrap = 8;      ///< Niveles de extrapolación (K_MAX)
        double safety   = 0.25;     ///< Factor de seguridad para ajuste de paso
        double initial_dtau = 1e-4; ///< Paso ficticio inicial
        double min_dtau = 1e-12;    ///< Paso ficticio minimo permitido
        double max_dtau = 1e-2;     ///< Paso ficticio maximo (conservador para KS)
        double energy_tol = 0.05;   ///< Tolerancia de energía: rechaza paso si |ΔH/H₀| > energy_tol
        bool   verbose  = false;    ///< Salida de diagnóstico
    };

    explicit Chain3BSIntegrator(const BSParameters& params = {});

    // -----------------------------------------------------------------------
    // INTERFAZ PRINCIPAL
    // -----------------------------------------------------------------------

    /**
     * @brief Integra con control de error estricto hasta t_target (tiempo físico).
     *
     * El control de paso se hace en tiempo ficticio τ.
     * cm_time se propaga a través de deriv.cm_time = g en modified_midpoint_step,
     * sin asignación manual, preservando la consistencia con la regularización.
     */
    void integrate_precise(
        Chain3State&       state,
        double             t_target,
        double&            t_achieved,
        const NBodySystem& full_system
    );

    /**
     * @brief Un paso BS: integra dtau_total en tiempo ficticio.
     *
     * @param state        Estado entrada/salida (solo modificado si accepted).
     * @param dtau_total   Paso ficticio a intentar.
     * @param error_estimate Error estimado (salida, siempre).
     * @param dtau_out     Paso sugerido para el siguiente intento (salida).
     * @param full_system  Sistema completo.
     * @return true si el paso fue aceptado.
     */
    bool bs_step(
        Chain3State&       state,
        double             dtau_total,
        double&            error_estimate,
        double&            dtau_out,
        const NBodySystem& full_system
    );

    /// Calcula la energía física directamente desde las variables KS del estado
    double compute_energy_from_state(const Chain3State& state) const;

    /// Error relativo de energía respecto a una referencia
    double compute_energy_error(const Chain3State& state, double E_ref) const;

    /// Energía física estándar (0.5*T_mikkola + U) — coincide con NBodySystem::total_energy()
    double compute_classical_energy(const Chain3State& state) const;  // energía clásica (½·T_Mikkola + U), ≠ state.energy

private:
    BSParameters bs_params_;

    // Secuencia de subdivisiones n_k = 2, 4, 6, 8, 10, 12, 14, 16
    static constexpr int K_MAX = 8;
    static constexpr std::array<int, K_MAX> n_seq_ = {2, 4, 6, 8, 10, 12, 14, 16};

    // -----------------------------------------------------------------------
    // MÉTODOS INTERNOS
    // -----------------------------------------------------------------------

    /**
     * @brief Paso medio modificado de Gragg.
     *
     * Propaga cm_time a través de deriv.cm_time = g (NO manualmente).
     */
    void modified_midpoint_step(
        const Chain3State& y0,
        double             dtau,      ///< Paso ficticio de cada sub-paso
        int                n_steps,
        const NBodySystem& full_system,
        Chain3State&       result
    );

    /// Extrapolación de Richardson (Neville-Aitken) de orden k_current
    void polynomial_extrapolation(
        const std::vector<Chain3State>& results,
        const std::vector<double>&      step_sizes,
        int                             k_current,
        Chain3State&                    extrapolated
    ) const;

    /// Error relativo entre dos niveles de extrapolación
    double estimate_error_bs(
        const Chain3State& high_order,
        const Chain3State& low_order
    ) const;

    /// Factor de escala para el siguiente paso según error y orden
    double step_scale_factor(double error, int k_opt) const;
};