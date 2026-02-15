#pragma once
#include "integrator.h"
#include "nbody_system.h"

class LeapfrogIntegrator : public Integrator {
public:
    void step(
        NBodySystem& system,
        double dt,
        const std::vector<bool>& used
    ) override {
        const size_t N = system.bodies.size();

        // === kick (half) ===
        auto acc = system.compute_accelerations();
        for (size_t i = 0; i < N; ++i) {
            if (used[i]) continue;
            system.bodies[i].velocity += (dt * 0.5) * acc[i];
        }

        // === drift ===
        for (size_t i = 0; i < N; ++i) {
            if (used[i]) continue;
            system.bodies[i].position += dt * system.bodies[i].velocity;
        }

        // === kick (half) ===
        acc = system.compute_accelerations();
        for (size_t i = 0; i < N; ++i) {
            if (used[i]) continue;
            system.bodies[i].velocity += (dt * 0.5) * acc[i];
        }
    }
};

// Nota: Esto no es un Leapfrog “perfecto” cuando hay exclusiones, pero:
// no integra 2 veces
// mantiene coherencia temporal
// permite seguir desarrollando el hibrido

