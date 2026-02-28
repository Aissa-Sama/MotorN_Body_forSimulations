// tests/test_chain_ks.cpp
#include <iostream>
#include <iomanip>
#include "vec3.h"
#include "vec4.h"
#include "chain_ks.h"

int main() {
    std::cout << std::setprecision(15);
    std::cout << "=== Test de transformaciones KS ===\n";

    // Probar con un vector aleatorio
    Vec3 r(1.2, 3.4, 5.6);  // Ahora funciona por el constructor
    std::cout << "r original: " << r.x << ", " << r.y << ", " << r.z << "\n";

    Vec4 u = r_to_ks(r);
    std::cout << "u = " << u.x << ", " << u.y << ", " << u.z << ", " << u.w << "\n";
    std::cout << "|u|^2 = " << u.norm2() << " (debe ser ~ |r| = " << r.norm() << ")\n";

    Vec3 r2 = ks_to_r(u);
    std::cout << "r recuperado: " << r2.x << ", " << r2.y << ", " << r2.z << "\n";
    std::cout << "Error posición: " << (r2 - r).norm() << "\n";

    // Probar con velocidades
    Vec3 v(0.1, 0.2, 0.3);
    Vec4 w = v_to_ks(v, u);
    Vec3 v2 = ks_to_v(u, w);
    std::cout << "v original: " << v.x << ", " << v.y << ", " << v.z << "\n";
    std::cout << "v recuperada: " << v2.x << ", " << v2.y << ", " << v2.z << "\n";
    std::cout << "Error velocidad: " << (v2 - v).norm() << "\n";

    return 0;
}