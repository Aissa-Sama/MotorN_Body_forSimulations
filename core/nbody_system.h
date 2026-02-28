// core/nbody_system.h - DEFINICIÓN CORRECTA
#pragma once
#include <vector>
#include "vec3.h"
#include "body.h"

class NBodySystem {
public:
    std::vector<Body> bodies;
    double G = 1.0;
    
    // Cache de aceleraciones
    mutable std::vector<Vec3> accelerations;
    mutable bool accelerations_valid = false;
    
    std::vector<Vec3> compute_accelerations() const;
    void invalidate_accelerations();
    
    // Invariantes físicos
    double kinetic_energy() const;
    double potential_energy() const;
    double total_energy() const;
    Vec3 total_momentum() const;
    Vec3 total_angular_momentum() const;
};