#pragma once
#include "vec3.h"

// ============================================================================
// BODY_TYPE — clasificación física del cuerpo
//
// Determina qué criterios de encuentro cercano aplican:
//   POINT_MASS       → sin radio físico. Solo interacción gravitacional pura.
//   STAR             → radio estelar conocido. Criterios TDE + colisión física.
//   COMPACT_OBJECT   → BH o NS. r_Schwarzschild como límite inferior. TDE activo.
//   PLANET           → radio planetario. Colisión física. TDE despreciable.
//   SATELLITE        → radio pequeño. Colisión física con planeta/estrella.
//
// Default: POINT_MASS — backward compatible con todo el código existente.
// ============================================================================
enum class BodyType {
    POINT_MASS      = 0,
    STAR            = 1,
    COMPACT_OBJECT  = 2,   // BH, NS, WD
    PLANET          = 3,
    SATELLITE       = 4
};

struct Body {
    Vec3   position;
    Vec3   velocity;
    double mass;

    // ── Extensiones Fase 7C ──────────────────────────────────────────────────
    // radius: radio físico en N-body units.
    //   - STAR:           R_sun en N-body units (típico: ~5e-3 para unidades pc)
    //   - COMPACT_OBJECT: r_Schwarzschild = 2GM/c² (en N-body units)
    //   - PLANET:         radio planetario
    //   - POINT_MASS:     0.0 (sin radio físico — masa puntual pura)
    //
    // Si radius == 0.0 → el cuerpo es tratado como masa puntual para todos los
    // criterios de colisión. Compatible con todo el código previo a Fase 7C.
    double   radius   = 0.0;
    BodyType type     = BodyType::POINT_MASS;

    // ── Constructores de conveniencia ─────────────────────────────────────────

    // Constructor legacy (compatible con código existente)
    Body() : position(), velocity(), mass(0.0), radius(0.0), type(BodyType::POINT_MASS) {}

    // Constructor con inicialización aggregate (posición, velocidad, masa)
    // — mantiene compatibilidad con {{pos}, {vel}, mass} en tests existentes
    Body(Vec3 p, Vec3 v, double m,
         double r = 0.0, BodyType t = BodyType::POINT_MASS)
        : position(p), velocity(v), mass(m), radius(r), type(t) {}

    // ── Helpers ───────────────────────────────────────────────────────────────

    bool is_extended()      const { return radius > 0.0; }
    bool is_compact()       const { return type == BodyType::COMPACT_OBJECT; }
    bool is_star()          const { return type == BodyType::STAR; }
    bool is_point_mass()    const { return type == BodyType::POINT_MASS; }

    // Radio de Schwarzschild en las mismas unidades que radius.
    // Significativo solo para COMPACT_OBJECT. c debe pasarse en N-body units.
    // r_S = 2GM/c²  con G=1 → r_S = 2*mass/c²
    double schwarzschild_radius(double c) const {
        return 2.0 * mass / (c * c);
    }
};
