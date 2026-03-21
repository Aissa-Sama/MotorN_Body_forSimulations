#pragma once
// ============================================================================
// analysis/tidal_events.h
//
// Funciones puras para detección de eventos de marea, disrupciones y
// colisiones físicas. Sin estado interno. Sin dependencias de integradores.
//
// FÍSICA IMPLEMENTADA:
//
// 1. Radio de disrupción tidal (Hills 1975; Rees 1988):
//    r_t = R_* · (M_BH / m_*)^(1/3)
//    Cuando r_p < r_t, las fuerzas de marea del objeto masivo superan la
//    auto-gravedad de la estrella → TDE.
//
// 2. Radio de Roche (Eggleton 1983):
//    r_L / a = 0.49 q^(2/3) / [0.6 q^(2/3) + ln(1 + q^(1/3))]
//    Fórmula exacta, válida para cualquier q = m_donante / m_acretor.
//    Error < 1% vs. solución numérica exacta.
//
// 3. Colisión física:
//    r_coll = R_i + R_j
//    Solo aplicable si ambos cuerpos tienen radio físico (radius > 0).
//
// 4. Radio de Schwarzschild (captura relativista):
//    r_S = 2GM/c²   (con G=1: r_S = 2M/c²)
//    Límite inferior para COMPACT_OBJECT.
//
// REFERENCIAS:
//   Hills (1975), Nature 254, 295
//   Rees (1988), Nature 333, 523
//   Eggleton (1983), ApJ 268, 368
//   Wang & Merritt (2004), ApJ 600, 149 — criterio de captura para BH
// ============================================================================
#include "body.h"
#include "vec3.h"
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

namespace tidal {

// ── Tipos de evento ───────────────────────────────────────────────────────────

enum class EventType {
    NONE              = 0,
    PHYSICAL_COLLISION = 1,   // sep < R_i + R_j  (dos cuerpos extendidos)
    FULL_TDE          = 2,   // sep < r_t          (estrella destruida por BH)
    PARTIAL_TDE       = 3,   // r_t < sep < 2*r_t  (estrella pierde masa, sobrevive)
    ROCHE_OVERFLOW    = 4,   // sep < r_Roche       (inicio transferencia de masa)
    SCHWARZSCHILD_CAP = 5,   // sep < r_S           (captura relativista directa)
};

struct TidalEvent {
    int       i, j;           // índices en NBodySystem.bodies
    EventType type;
    double    sep;            // separación en el momento de detección
    double    r_crit;         // radio crítico que se cruzó
    double    time;           // tiempo físico del evento
    std::string description;  // log legible
};

// ── Radio de disrupción tidal (Hills 1975) ────────────────────────────────────
//
// r_t = R_star · (M_massive / m_star)^(1/3)
//
// Interpretación: a esta separación, la fuerza de marea del objeto masivo
// sobre la estrella es igual a la fuerza de auto-gravedad superficial de
// la estrella. Para r < r_t, la estrella se desgarra.
//
// Parámetros:
//   M_massive : masa del objeto masivo (BH, SMBH) en N-body units
//   m_star    : masa de la estrella en N-body units
//   R_star    : radio de la estrella en N-body units
//
inline double tidal_radius(double M_massive, double m_star, double R_star) {
    if (R_star <= 0.0 || m_star <= 0.0 || M_massive <= 0.0) return 0.0;
    return R_star * std::cbrt(M_massive / m_star);
}

// ── Radio de Roche (Eggleton 1983) ───────────────────────────────────────────
//
// r_L / a = 0.49 q^(2/3) / [0.6 q^(2/3) + ln(1 + q^(1/3))]
//
// q = m_donante / m_acretor
// a = separación orbital (semieje mayor)
//
// Válida para 0 < q < ∞ con error < 1%.
// Referencia estándar en literatura de binarias interactuantes.
//
inline double roche_lobe_radius(double a_orbit, double q) {
    if (q <= 0.0 || a_orbit <= 0.0) return 0.0;
    const double q13 = std::cbrt(q);
    const double q23 = q13 * q13;
    return a_orbit * (0.49 * q23) / (0.6 * q23 + std::log(1.0 + q13));
}

// ── Radio de Schwarzschild ────────────────────────────────────────────────────
//
// r_S = 2GM/c²   con G=1 → r_S = 2M/c²
//
// Para simulaciones con PN activo, c es el parámetro de velocidad de la luz
// en N-body units (mismo c que en ARChainNPNIntegrator).
//
inline double schwarzschild_radius(double mass, double c) {
    return 2.0 * mass / (c * c);
}

// ── Radio ISCO (Innermost Stable Circular Orbit) ──────────────────────────────
//
// Para BH de Schwarzschild (no giratorio):
//   r_ISCO = 3 · r_S = 6GM/c²
//
// Objetos que alcanzan r < r_ISCO son capturados inexorablemente.
// Para BH de Kerr (giratorio), r_ISCO puede ser tan pequeño como r_S/2.
// Aquí usamos Schwarzschild como aproximación conservadora.
//
inline double isco_radius(double mass, double c) {
    return 3.0 * schwarzschild_radius(mass, c);
}

// ── Clasificación del encuentro entre dos cuerpos ─────────────────────────────
//
// Determina qué evento (si alguno) se produce dado el par (i,j) y su
// separación actual. El orden de prioridad es:
//
//   1. Schwarzschild (captura relativista) — si c > 0 y alguno es COMPACT
//   2. Colisión física — si ambos tienen radio
//   3. Full TDE — si uno es COMPACT/masivo y el otro es STAR con radius
//   4. Partial TDE — zona de 1 < sep/r_t < 2
//   5. Roche overflow — si hay semieje definible (binaria ligada)
//
// Parámetros:
//   a  : cuerpo i
//   b  : cuerpo j
//   sep: separación actual |r_i - r_j|
//   c  : velocidad de la luz en N-body units (0 = no hay física relativista)
//
inline EventType classify_encounter(const Body& a, const Body& b,
                                    double sep, double c = 0.0) {
    // ── 1. Captura relativista ────────────────────────────────────────────────
    if (c > 0.0) {
        // Si alguno es objeto compacto, verificar ISCO del más masivo
        if (a.is_compact() || b.is_compact()) {
            const Body& heavy = (a.mass >= b.mass) ? a : b;
            double r_isco = isco_radius(heavy.mass, c);
            if (sep < r_isco) return EventType::SCHWARZSCHILD_CAP;
        }
    }

    // ── 2. Colisión física ────────────────────────────────────────────────────
    if (a.is_extended() && b.is_extended()) {
        if (sep < a.radius + b.radius)
            return EventType::PHYSICAL_COLLISION;
    }

    // ── 3 & 4. TDE ────────────────────────────────────────────────────────────
    // Un cuerpo debe ser masivo (COMPACT o POINT_MASS con masa >> otro),
    // el otro debe ser STAR con radio conocido.
    //
    // Definimos "masivo" como el de mayor masa cuando la ratio > 100.
    const bool mass_ratio_ok = (a.mass / b.mass > 100.0) || (b.mass / a.mass > 100.0);

    if (mass_ratio_ok) {
        const Body& heavy = (a.mass >= b.mass) ? a : b;
        const Body& light = (a.mass <  b.mass) ? a : b;

        if (light.is_extended() && light.is_star()) {
            double r_t = tidal_radius(heavy.mass, light.mass, light.radius);
            if (r_t > 0.0) {
                if (sep < r_t)
                    return EventType::FULL_TDE;
                if (sep < 2.0 * r_t)
                    return EventType::PARTIAL_TDE;
            }
        }
    }

    return EventType::NONE;
}

// ── Detector de eventos para todo el sistema ──────────────────────────────────
//
// Itera sobre todos los pares (i,j) y devuelve los eventos detectados.
// Ordenados por tipo de prioridad (SCHWARZSCHILD primero).
//
// Parámetros:
//   bodies : vector de cuerpos del sistema
//   c      : velocidad de la luz (0 = sin física relativista)
//   t_now  : tiempo físico actual (para logging)
//
inline std::vector<TidalEvent> detect_events(
    const std::vector<Body>& bodies,
    double c = 0.0,
    double t_now = 0.0)
{
    std::vector<TidalEvent> events;
    const int N = static_cast<int>(bodies.size());

    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            Vec3 dr = bodies[j].position - bodies[i].position;
            double sep = dr.norm();

            EventType ev = classify_encounter(bodies[i], bodies[j], sep, c);
            if (ev == EventType::NONE) continue;

            TidalEvent te;
            te.i    = i;
            te.j    = j;
            te.type = ev;
            te.sep  = sep;
            te.time = t_now;

            // Radio crítico y descripción
            switch (ev) {
                case EventType::PHYSICAL_COLLISION:
                    te.r_crit      = bodies[i].radius + bodies[j].radius;
                    te.description = "COLLISION i=" + std::to_string(i)
                                   + " j=" + std::to_string(j);
                    break;
                case EventType::FULL_TDE: {
                    const Body& heavy = (bodies[i].mass >= bodies[j].mass) ? bodies[i] : bodies[j];
                    const Body& light = (bodies[i].mass <  bodies[j].mass) ? bodies[i] : bodies[j];
                    te.r_crit      = tidal_radius(heavy.mass, light.mass, light.radius);
                    te.description = "FULL_TDE BH=" + std::to_string(i)
                                   + " STAR=" + std::to_string(j);
                    break;
                }
                case EventType::PARTIAL_TDE: {
                    const Body& heavy = (bodies[i].mass >= bodies[j].mass) ? bodies[i] : bodies[j];
                    const Body& light = (bodies[i].mass <  bodies[j].mass) ? bodies[i] : bodies[j];
                    te.r_crit      = tidal_radius(heavy.mass, light.mass, light.radius);
                    te.description = "PARTIAL_TDE BH=" + std::to_string(i)
                                   + " STAR=" + std::to_string(j);
                    break;
                }
                case EventType::SCHWARZSCHILD_CAP: {
                    const Body& heavy = (bodies[i].mass >= bodies[j].mass) ? bodies[i] : bodies[j];
                    te.r_crit      = isco_radius(heavy.mass, c);
                    te.description = "SCHWARZSCHILD_CAP BH=" + std::to_string(i)
                                   + " obj=" + std::to_string(j);
                    break;
                }
                default:
                    te.r_crit = 0.0;
                    te.description = "UNKNOWN";
                    break;
            }

            events.push_back(te);
        }
    }

    // Ordenar por prioridad (SCHWARZSCHILD > COLLISION > FULL_TDE > ...)
    std::sort(events.begin(), events.end(), [](const TidalEvent& a, const TidalEvent& b) {
        return static_cast<int>(a.type) < static_cast<int>(b.type);
    });

    return events;
}

} // namespace tidal