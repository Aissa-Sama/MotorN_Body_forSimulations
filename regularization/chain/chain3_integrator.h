// regularization/chain/chain3_integrator.h
#pragma once
#include "chain_state.h"
#include "nbody_system.h"
#include <memory>

/**
 * @brief Estructura para parámetros de integración
 */
struct IntegrationParams {
    double abs_tol = 1e-12;      // Tolerancia absoluta para error
    double rel_tol = 1e-12;      // Tolerancia relativa para error
    double min_dtau = 1e-8;       // Paso ficticio mínimo
    double max_dtau = 1e-1;       // Paso ficticio máximo
    double safety_factor = 0.9;   // Factor de seguridad para ajuste de paso
    bool verbose = false;         // Modo verbose para debugging
};

/**
 * @brief Integrador de cadena para sistemas de 3 cuerpos.
 * 
 * Implementa la regularización en cadena según Mikkola & Aarseth (1993).
 * Esta clase se encarga de:
 * - Transformar entre coordenadas físicas y variables de cadena.
 * - Integrar el sistema en tiempo ficticio.
 * - Reconstruir el sistema físico.
 */
class Chain3Integrator {
public:
    /**
     * @brief Constructor.
     * @param internal_dt Paso de tiempo ficticio interno (valor típico: 1e-4).
     */
    explicit Chain3Integrator(double internal_dt);

    /**
     * @brief Inicializa el estado de la cadena a partir del sistema físico.
     * 
     * @param system Sistema N-body completo.
     * @param idx1 Índice del primer cuerpo en la cadena.
     * @param idx2 Índice del segundo cuerpo en la cadena.
     * @param idx3 Índice del tercer cuerpo en la cadena.
     * @return Chain3State Estado inicial de la cadena.
     * 
     * @pre Los tres cuerpos deben existir en system.
     * @post El estado de la cadena es consistente con las posiciones y velocidades.
     */
    Chain3State initialize(const NBodySystem& system, int idx1, int idx2, int idx3);

    /**
     * @brief Reconstruye el sistema físico a partir del estado de la cadena.
     * 
     * @param state Estado actual de la cadena.
     * @param system Sistema N-body (se modificarán los cuerpos idx1, idx2, idx3).
     * @param idx1 Índice del primer cuerpo.
     * @param idx2 Índice del segundo cuerpo.
     * @param idx3 Índice del tercer cuerpo.
     */
    void write_back(const Chain3State& state, NBodySystem& system, 
                    int idx1, int idx2, int idx3);

    /**
     * @brief Versión simple de integración (sin control de error).
     * @param state Estado de la cadena (entrada/salida).
     * @param dt_global Paso de tiempo físico a integrar.
     * @param full_system Sistema completo (para perturbaciones futuras).
     */
    void integrate_simple(Chain3State& state, double dt_global, 
                          const NBodySystem& full_system);

    /**
     * @brief Integra la cadena desde el tiempo actual hasta t_target (tiempo físico).
     * 
     * @param state Estado de la cadena (entrada/salida).
     * @param t_target Tiempo físico objetivo a alcanzar.
     * @param t_achieved Tiempo físico realmente avanzado (salida).
     * @param params Parámetros de integración.
     * @param full_system Sistema completo (para perturbaciones futuras).
     */
    void integrate(Chain3State& state,
                   double t_target,
                   double& t_achieved,
                   const IntegrationParams& params,
                   const NBodySystem& full_system);

    // ============================================================================
    // MÉTODOS PÚBLICOS PARA TESTS (ANTES ERAN PRIVADOS)
    // ============================================================================
    
    /**
     * @brief Calcula las derivadas de las variables KS (público para tests).
     * @param state Estado actual.
     * @param full_system Sistema completo (para perturbaciones futuras).
     * @param deriv Estructura donde se guardarán las derivadas.
     */
    void compute_derivatives(const Chain3State& state, 
                             const NBodySystem& full_system,
                             Chain3State& deriv) const;
    
    /**
     * @brief Realiza un paso RK4 completo (público para tests).
     * @param state Estado de entrada/salida.
     * @param dtau Paso de tiempo ficticio.
     * @param full_system Sistema completo.
     */
    void rk4_step(Chain3State& state, double dtau, const NBodySystem& full_system);
    
    /**
     * @brief Actualiza el tiempo físico (público para tests).
     * @param state Estado actual.
     * @param dtau Paso de tiempo ficticio.
     */
    void update_time(Chain3State& state, double dtau);

private:
    double internal_dt;  ///< Paso de tiempo ficticio interno

    // Estos métodos permanecen privados (solo uso interno)
    void rk4_step_with_error(Chain3State& state, double dtau, 
                             const NBodySystem& full_system,
                             Chain3State& error_estimate);
    
    double estimate_error(const Chain3State& state1, const Chain3State& state2);
};