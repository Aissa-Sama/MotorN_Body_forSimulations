// regularization/chain/chain_ks.h
#pragma once
#include "vec4.h"
#include "vec3.h"
#include <stdexcept>
#include <cmath>

// ============================================================================
// MATRIZ DE LEVI-CIVITA Y OPERACIONES
// ============================================================================

inline void levi_civita_matrix(const Vec4& u, double L[3][4]) {
    L[0][0] =  u.x;  L[0][1] = -u.y;  L[0][2] = -u.z;  L[0][3] =  u.w;
    L[1][0] =  u.y;  L[1][1] =  u.x;  L[1][2] = -u.w;  L[1][3] = -u.z;
    L[2][0] =  u.z;  L[2][1] =  u.w;  L[2][2] =  u.x;  L[2][3] =  u.y;
}

/**
 * @brief Multiplica L(u) por vector v (4D) → resultado 3D
 */
inline Vec3 levi_civita_multiply(const Vec4& u, const Vec4& v) {
    return Vec3{
        u.x * v.x - u.y * v.y - u.z * v.z + u.w * v.w,
        u.y * v.x + u.x * v.y - u.w * v.z - u.z * v.w,
        u.z * v.x + u.w * v.y + u.x * v.z + u.y * v.w
    };
}

/**
 * @brief Multiplica L(u)^T por vector r (3D) → resultado 4D
 */
inline Vec4 levi_civita_transpose_multiply(const Vec4& u, const Vec3& r) {
    return Vec4{
        u.x * r.x + u.y * r.y + u.z * r.z,
        -u.y * r.x + u.x * r.y + u.w * r.z,
        -u.z * r.x - u.w * r.y + u.x * r.z,
        u.w * r.x - u.z * r.y + u.y * r.z
    };
}

// ============================================================================
// TRANSFORMACIONES KS
// ============================================================================

inline Vec3 ks_to_r(const Vec4& u) {
    return levi_civita_multiply(u, u);
}

/**
 * @brief Transformación KS inversa completa y robusta
 */
inline Vec4 r_to_ks(const Vec3& r) {
    double R = r.norm();
    if (R < 1e-12) throw std::runtime_error("r_to_ks: |r| too small");
    
    if (r.x >= 0.0) {
        double t = std::sqrt((R + r.x) / 2.0);
        return Vec4(t, r.y/(2.0*t), r.z/(2.0*t), 0.0);
    } else {
        double t = std::sqrt((R - r.x) / 2.0);
        return Vec4(r.y/(2.0*t), t, 0.0, r.z/(2.0*t));
    }
}

// ============================================================================
// TRANSFORMACIONES DE MOMENTOS
// ============================================================================

inline Vec4 w_to_ks(const Vec3& W, const Vec4& u) {
    return levi_civita_transpose_multiply(u, W) * 2.0;
}

inline Vec3 ks_to_w(const Vec4& u, const Vec4& P) {
    double u2 = u.norm2();
    if (u2 < 1e-24) throw std::runtime_error("ks_to_w: |u|^2 too small");
    return levi_civita_multiply(u, P) * (0.5 / u2);
}

// Alias para compatibilidad
inline Vec3 ks_to_v(const Vec4& u, const Vec4& p) { return ks_to_w(u, p); }
inline Vec4 v_to_ks(const Vec3& v, const Vec4& u) { return w_to_ks(v, u); }

/**
 * @brief Multiplica L(u) por vector v (3D tratado como 4D con w=0) → resultado 3D
 * 
 * Útil para calcular términos como L(P) * A donde A es 3D
 */
inline Vec3 levi_civita_multiply_vec3(const Vec4& u, const Vec3& v) {
    // Tratar v como (v.x, v.y, v.z, 0)
    return Vec3{
        u.x * v.x - u.y * v.y - u.z * v.z,  // + u.w * 0
        u.y * v.x + u.x * v.y - u.z * 0   - u.w * v.z,
        u.z * v.x + u.w * v.y + u.x * v.z + u.y * 0
    };
}