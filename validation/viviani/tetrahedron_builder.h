#pragma once
#include "vec3.h"
#include "nbody_system.h"

class VivianiTetrahedron {
public:
    struct Tetrahedron {
        Vec3 vertices[4];        // V1, V2, V3, V4
        Vec3 face_centroids[4];  // B1=P1, B2=P2, B3=P3, B4
        Vec3 centroid;           // Gc = (B1+B2+B3+B4)/4
        double volume;
        double alpha;
    };

    // B4 = G + α·L̂ (ecuación del paper §4.3)
    static Tetrahedron build(const NBodySystem& system, int idx1, int idx2, int idx3, 
                            double epsilon = 0.1);
    
    // Verificar no degeneración: det(B2-B1, B3-B1, B4-B1) ≠ 0
    static bool check_non_degeneracy(const Tetrahedron& tet, double tol = 1e-12);
};