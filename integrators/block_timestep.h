// integrators/block_timestep.h
// FASE 7B — Block timestep (Aarseth 2003, cap. 2)
//
// PRINCIPIO:
//   Cada cuerpo libre (LEAF) avanza con su propio paso de tiempo basado en
//   su escala de tiempo dinamica local. Los pasos son potencias de 2 del
//   paso minimo para garantizar sincronizacion exacta entre cuerpos.
//
// CRITERIO DE PASO (Aarseth 2003, ec. 2.5):
//   dt_i = eta * sqrt(|a_i| / |jerk_i|)
//
//   donde:
//     a_i    = aceleracion gravitacional del cuerpo i
//     jerk_i = da_i/dt = derivada temporal de la aceleracion
//     eta    = parametro de precision (tipico: 0.02 - 0.04)
//
// BLOCK TIMESTEP:
//   El paso continuo se redondea a la potencia de 2 mas cercana hacia abajo:
//     k_i    = floor(log2(dt_i / dt_min))
//     k_i    = clamp(k_i, 0, k_max)
//     dt_i   = dt_min * 2^k_i
//
//   Esto garantiza que cualquier par de cuerpos se sincronice exactamente
//   cada dt_max = dt_min * 2^k_max pasos.
//
// JERK GRAVITACIONAL:
//   jerk_i = Σ_{j≠i} G m_j [ v_ij/r_ij³ - 3(v_ij·r_ij) r_ij / r_ij^5 ]
//
//   Identico a la formula en gw_observables.h (compute_quadrupole_dddot).
//   Se recalcula aqui para no crear dependencia circular con analysis/.
//
// INTEGRACION EN LA ARQUITECTURA:
//   HierarchicalIntegrator::step() pasa dt_max al campo lejano.
//   BlockTimestep::assign() devuelve un vector<double> con el dt asignado
//   a cada cuerpo. Los LEAFs usan su dt individual; los subsistemas
//   regularizados (PAIR_KS, GROUP_AR_CHAIN) usan dt_max sin cambio
//   porque ya tienen control de paso interno.
//
// REFERENCIAS:
//   Aarseth, S.J. (2003). Gravitational N-Body Simulations. Cambridge.
//   Makino, J. & Aarseth, S.J. (1992). PASJ 44, 141 — criterio de paso.
// ============================================================================
#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <algorithm>
#include "nbody_system.h"
#include "vec3.h"

class BlockTimestep {
public:

    struct Params {
        double eta    = 0.03;   ///< Factor de precision. Mayor = mas preciso pero mas lento.
        double dt_min = 1e-5;   ///< Paso minimo: nunca menor que esto.
        double dt_max = 0.1;    ///< Paso maximo: nunca mayor que esto (= dt del campo lejano).
        int    k_max  = 10;     ///< Numero maximo de niveles de bloque (2^k_max * dt_min <= dt_max).
    };

    explicit BlockTimestep(const Params& p = Params{}) : params_(p) {}

    // ── Interfaz principal ───────────────────────────────────────────────────

    /// Calcula el jerk gravitacional de todos los cuerpos.
    /// jerk[i] = da_i/dt = Σ_{j≠i} G m_j [v_ij/r³ - 3(v_ij·r) r/r^5]
    static std::vector<Vec3> compute_jerks(
        const NBodySystem& sys,
        const std::vector<Vec3>& accs);

    /// Calcula el dt de bloque para un cuerpo dado su |a| y |jerk|.
    /// Retorna dt_min * 2^k con k en [0, k_max].
    double block_dt(double acc_norm, double jerk_norm) const;

    /// Asigna un dt de bloque a cada cuerpo del sistema.
    /// used[i] = true → cuerpo en subsistema regularizado → usa dt_max.
    /// used[i] = false → LEAF → usa su dt individual.
    ///
    /// El dt_max se usa como techo (ningun cuerpo avanza mas que dt_max).
    std::vector<double> assign(
        const NBodySystem& sys,
        const std::vector<Vec3>& accs,
        const std::vector<Vec3>& jerks,
        const std::vector<bool>& used,
        double dt_max) const;

    const Params& params() const { return params_; }
    void set_params(const Params& p) { params_ = p; }

private:
    Params params_;
};

// ============================================================================
// IMPLEMENTACION INLINE
// ============================================================================

inline std::vector<Vec3> BlockTimestep::compute_jerks(
    const NBodySystem& sys,
    const std::vector<Vec3>& accs)
{
    const int N = static_cast<int>(sys.bodies.size());
    std::vector<Vec3> jerks(N);

    for (int i = 0; i < N; ++i) {
        Vec3 jk;
        for (int j = 0; j < N; ++j) {
            if (i == j) continue;
            Vec3 rij = sys.bodies[j].position - sys.bodies[i].position;
            Vec3 vij = sys.bodies[j].velocity - sys.bodies[i].velocity;
            const double r2 = rij.norm2();
            const double r  = std::sqrt(r2);
            const double r3 = r2 * r;
            const double r5 = r3 * r2;
            if (r3 < 1e-90) continue;
            const double v_dot_r = dot(vij, rij);
            const double fac = sys.G * sys.bodies[j].mass;
            jk = jk + vij * (fac / r3) - rij * (3.0 * fac * v_dot_r / r5);
        }
        jerks[i] = jk;
    }
    return jerks;
}

inline double BlockTimestep::block_dt(double acc_norm, double jerk_norm) const
{
    // Criterio de Aarseth: dt = eta * sqrt(|a| / |jerk|)
    if (jerk_norm < 1e-30 || acc_norm < 1e-30)
        return params_.dt_max;

    const double dt_crit = params_.eta * std::sqrt(acc_norm / jerk_norm);

    // Redondear a potencia de 2 hacia abajo
    const double ratio = dt_crit / params_.dt_min;
    if (ratio < 1.0) return params_.dt_min;

    const int k = std::min(
        static_cast<int>(std::floor(std::log2(ratio))),
        params_.k_max);

    return params_.dt_min * std::pow(2.0, static_cast<double>(k));
}

inline std::vector<double> BlockTimestep::assign(
    const NBodySystem& sys,
    const std::vector<Vec3>& accs,
    const std::vector<Vec3>& jerks,
    const std::vector<bool>& used,
    double dt_max) const
{
    const int N = static_cast<int>(sys.bodies.size());
    std::vector<double> dts(N, dt_max);

    for (int i = 0; i < N; ++i) {
        if (used[i]) {
            // Cuerpo en subsistema regularizado — usa dt_max
            dts[i] = dt_max;
            continue;
        }
        const double an = accs[i].norm();
        const double jn = jerks[i].norm();
        // Aplicar dt_max como techo
        dts[i] = std::min(block_dt(an, jn), dt_max);
    }
    return dts;
}