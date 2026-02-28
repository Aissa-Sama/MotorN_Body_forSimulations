// core/vec3.h
#pragma once
#include <cmath>
#include <iostream>

struct Vec3 {
    double x, y, z;

    // Constructores
    Vec3() : x(0.0), y(0.0), z(0.0) {}
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    // Operador unario negativo
    Vec3 operator-() const {
        return { -x, -y, -z };
    }

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
        double inv = 1.0 / s;
        x *= inv; y *= inv; z *= inv;
        return *this;
    }

    // Norma
    double norm() const {
        return std::sqrt(x*x + y*y + z*z);
    }

    double norm2() const {
        return x*x + y*y + z*z;
    }

    // Producto punto
    double dot(const Vec3& other) const {
        return x*other.x + y*other.y + z*other.z;
    }

    // Imprimir
    void print(const std::string& name = "") const {
        if (!name.empty()) std::cout << name << " ";
        std::cout << "(" << x << ", " << y << ", " << z << ")\n";
    }
};

// Operadores libres
inline Vec3 operator*(double s, const Vec3& v) {
    return { s * v.x, s * v.y, s * v.z };
}

// Funciones libres para compatibilidad
inline double dot(const Vec3& a, const Vec3& b) {
    return a.dot(b);
}

inline double norm(const Vec3& v) {
    return v.norm();
}

// Producto vectorial (cross product)
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return Vec3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    );
}
