#pragma once
#include <vector>
#include "body.h"

class NBodySystem {
public:
    std::vector<Body> bodies;
    double G = 1.0;

    // Caché de aceleraciones
    mutable std::vector<Vec3> accelerations;  // mutable permite modificarlo en const
    mutable bool accelerations_valid = false;

    // =========================
    // Dinámica
    // =========================

    std::vector<Vec3> compute_accelerations() const {
        if (accelerations_valid) {
            return accelerations;
        }
        
        accelerations.assign(bodies.size(), {0, 0, 0});
        
        for (size_t i = 0; i < bodies.size(); ++i) {
            for (size_t j = 0; j < bodies.size(); ++j) {
                if (i == j) continue;
                
                Vec3 r = bodies[j].position - bodies[i].position;
                double d = norm(r) + 1e-9;
                
                accelerations[i] = accelerations[i] + 
                    (G * bodies[j].mass / (d * d * d)) * r;
            }
        }
        
        accelerations_valid = true;
        return accelerations;
    }
    
    void invalidate_accelerations() {
        accelerations_valid = false;
    }
    
    // =========================
    // Invariantes físicos
    // =========================
double kinetic_energy() const {
    double E = 0.0;
    for (const auto& b : bodies) {
        E += 0.5 * b.mass * dot(b.velocity, b.velocity);
    }
    return E;
}

double potential_energy() const {
    double E = 0.0;
    for (size_t i = 0; i < bodies.size(); ++i) {
        for (size_t j = i + 1; j < bodies.size(); ++j) {
            Vec3 r = bodies[j].position - bodies[i].position;
            double d = norm(r) + 1e-9;
            E -= G * bodies[i].mass * bodies[j].mass / d;
        }
    }
    return E;
}

double total_energy() const {
    return kinetic_energy() + potential_energy();
}

    Vec3 total_angular_momentum() const {
        Vec3 L{0, 0, 0};
        for (const auto& b : bodies) {
            // L = r × p = r × (m v)
            Vec3 p = b.mass * b.velocity;
            L.x += b.position.y * p.z - b.position.z * p.y;
            L.y += b.position.z * p.x - b.position.x * p.z;
            L.z += b.position.x * p.y - b.position.y * p.x;
        }
        return L;
    }
    
    Vec3 total_momentum() const {
        Vec3 P{0, 0, 0};
        for (const auto& b : bodies) {
            P = P + b.mass * b.velocity;
        }
        return P;
    }

};

// Nota: La implementación de potential_energy() está bien porque usa j = i+1, evitando la doble contabilidad.