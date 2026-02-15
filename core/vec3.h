#pragma once
#include <cmath>

struct Vec3 {
    double x, y, z;

    // Operadores miembro
    Vec3 operator+(const Vec3& other) const {
        return { x + other.x, y + other.y, z + other.z };
    }

    Vec3 operator-(const Vec3& other) const {
        return { x - other.x, y - other.y, z - other.z };
    }

    Vec3 operator*(double s) const {
        return { x * s, y * s, z * s };
    }

    Vec3 operator/(double s) const {
        return { x / s, y / s, z / s };
    }

    Vec3& operator+=(const Vec3& other) {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }

    Vec3& operator-=(const Vec3& other) {
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }

    Vec3& operator*=(double s) {
        x *= s; y *= s; z *= s;
        return *this;
    }

    Vec3& operator/=(double s) {
        x /= s; y /= s; z /= s;
        return *this;
    }
};

// Operadores libres
inline Vec3 operator*(double s, const Vec3& v) {
    return { s * v.x, s * v.y, s * v.z };
}

inline double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double norm(const Vec3& v) {
    return std::sqrt(dot(v, v));
}

// Con esto: Vec3 / double funciona
// Y con esto: double * Vec3 funciona; +=, *= etc. funcionan
// ks y binarios funcionan también.
