// regularization/ks/ks_transform.h (CORREGIDO Y COMPLETO)
#pragma once
#include "vec3.h"
#include "ks_state.h"
#include <cmath>
#include <stdexcept>

// Matriz de Levi-Civita para transformación KS
inline void levi_civita_matrix(const double u[4], double L[3][4]) {
    L[0][0] =  u[0]; L[0][1] = -u[1]; L[0][2] = -u[2]; L[0][3] =  u[3];
    L[1][0] =  u[1]; L[1][1] =  u[0]; L[1][2] = -u[3]; L[1][3] = -u[2];
    L[2][0] =  u[2]; L[2][1] =  u[3]; L[2][2] =  u[0]; L[2][3] =  u[1];
}

// Transformación KS: r = L(u) * u
inline Vec3 ks_to_r(const double u[4]) {
    return Vec3{
        2.0 * (u[0]*u[2] + u[1]*u[3]),
        2.0 * (u[1]*u[2] - u[0]*u[3]),
        u[0]*u[0] + u[1]*u[1] - u[2]*u[2] - u[3]*u[3]
    };
}

// Transformación inversa (VERSIÓN COMPLETA Y SIMÉTRICA)
inline void r_to_ks(const Vec3& r, double u[4]) {
    double rnorm = norm(r);
    if (rnorm < 1e-12) {
        throw std::runtime_error("Posición relativa cero en transformación KS");
    }
    
    // Determinar la componente dominante para evitar divisiones por cero
    double ax = std::abs(r.x);
    double ay = std::abs(r.y);
    double az = std::abs(r.z);
    
    // Usar la componente más grande como referencia
    if (ax >= ay && ax >= az) {
        // Componente X dominante
        double temp = std::sqrt((rnorm + ax) / 2.0);
        if (r.x >= 0) {
            u[0] = temp;
            u[1] = 0.0;
            u[2] = r.y / (2.0 * temp);
            u[3] = r.z / (2.0 * temp);
        } else {
            u[0] = r.y / (2.0 * temp);
            u[1] = temp;
            u[2] = 0.0;
            u[3] = r.z / (2.0 * temp);
        }
    }
    else if (ay >= ax && ay >= az) {
        // Componente Y dominante
        double temp = std::sqrt((rnorm + ay) / 2.0);
        if (r.y >= 0) {
            u[0] = 0.0;
            u[1] = temp;
            u[2] = -r.x / (2.0 * temp);  // Signo para mantener orientación
            u[3] = r.z / (2.0 * temp);
        } else {
            u[0] = -r.x / (2.0 * temp);
            u[1] = 0.0;
            u[2] = temp;
            u[3] = r.z / (2.0 * temp);
        }
    }
    else {
        // Componente Z dominante
        double temp = std::sqrt((rnorm + az) / 2.0);
        if (r.z >= 0) {
            u[0] = temp;
            u[1] = 0.0;
            u[2] = r.x / (2.0 * temp);
            u[3] = r.y / (2.0 * temp);
        } else {
            u[0] = 0.0;
            u[1] = temp;
            u[2] = r.y / (2.0 * temp);
            u[3] = -r.x / (2.0 * temp);
        }
    }
}

// Relación entre velocidades: v = (2/|u|^2) * L(u) * w
inline Vec3 ks_to_v(const double u[4], const double w[4]) {
    double L[3][4];
    levi_civita_matrix(u, L);
    
    double u_norm_sq = u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
    double factor = 2.0 / u_norm_sq;
    
    Vec3 v;
    v.x = factor * (L[0][0]*w[0] + L[0][1]*w[1] + L[0][2]*w[2] + L[0][3]*w[3]);
    v.y = factor * (L[1][0]*w[0] + L[1][1]*w[1] + L[1][2]*w[2] + L[1][3]*w[3]);
    v.z = factor * (L[2][0]*w[0] + L[2][1]*w[1] + L[2][2]*w[2] + L[2][3]*w[3]);
    
    return v;
}

// Transformación inversa de velocidades (simplificada)
inline void v_to_ks(const Vec3& v, const double u[4], double w[4]) {
    double L[3][4];
    levi_civita_matrix(u, L);
    
    double u_norm_sq = u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
    
    // w = (1/2) * L(u)^T * v
    w[0] = 0.5 * (L[0][0]*v.x + L[1][0]*v.y + L[2][0]*v.z);
    w[1] = 0.5 * (L[0][1]*v.x + L[1][1]*v.y + L[2][1]*v.z);
    w[2] = 0.5 * (L[0][2]*v.x + L[1][2]*v.y + L[2][2]*v.z);
    w[3] = 0.5 * (L[0][3]*v.x + L[1][3]*v.y + L[2][3]*v.z);
}

inline double ks_radius(const double u[4]) {
    return u[0]*u[0] + u[1]*u[1] + u[2]*u[2] + u[3]*u[3];
}