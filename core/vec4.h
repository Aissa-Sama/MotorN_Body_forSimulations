// core/vec4.h
#pragma once
#include <cmath>
#include <iostream>

/**
 * @brief Vector de 4 componentes para coordenadas KS.
 * 
 * Proporciona operadores aritméticos básicos y funciones
 * para calcular norma y producto punto.
 */
struct Vec4 {
    double x, y, z, w;

    // Constructores
    Vec4() : x(0.0), y(0.0), z(0.0), w(0.0) {}
    Vec4(double x_, double y_, double z_, double w_) 
        : x(x_), y(y_), z(z_), w(w_) {}

    // Operador unario negativo (¡NUEVO!)
    Vec4 operator-() const {
        return { -x, -y, -z, -w };
    }

    // Operadores aritméticos
    Vec4 operator+(const Vec4& other) const {
        return {x + other.x, y + other.y, z + other.z, w + other.w};
    }

    Vec4 operator-(const Vec4& other) const {
        return {x - other.x, y - other.y, z - other.z, w - other.w};
    }

    Vec4 operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar, w * scalar};
    }

    Vec4 operator/(double scalar) const {
        double inv = 1.0 / scalar;
        return {x * inv, y * inv, z * inv, w * inv};
    }

    // Operadores de asignación compuesta
    Vec4& operator+=(const Vec4& other) {
        x += other.x; y += other.y; z += other.z; w += other.w;
        return *this;
    }

    Vec4& operator-=(const Vec4& other) {
        x -= other.x; y -= other.y; z -= other.z; w -= other.w;
        return *this;
    }

    Vec4& operator*=(double scalar) {
        x *= scalar; y *= scalar; z *= scalar; w *= scalar;
        return *this;
    }

    Vec4& operator/=(double scalar) {
        double inv = 1.0 / scalar;
        x *= inv; y *= inv; z *= inv; w *= inv;
        return *this;
    }

    // Norma al cuadrado (útil para evitar sqrt innecesario)
    double norm2() const {
        return x*x + y*y + z*z + w*w;
    }

    // Norma euclidiana
    double norm() const {
        return std::sqrt(norm2());
    }

    // Producto punto
    double dot(const Vec4& other) const {
        return x*other.x + y*other.y + z*other.z + w*other.w;
    }

    // Imprimir (para debugging)
    void print(const std::string& name = "") const {
        if (!name.empty()) std::cout << name << " ";
        std::cout << "(" << x << ", " << y << ", " << z << ", " << w << ")\n";
    }
};

// Operadores libres para multiplicación escalar por vector
inline Vec4 operator*(double scalar, const Vec4& v) {
    return {v.x * scalar, v.y * scalar, v.z * scalar, v.w * scalar};
}