// integrators/euler_integrator.h
#pragma once
#include "integrator.h"

class EulerIntegrator : public Integrator {
public:
    void step(NBodySystem& system, double dt, const std::vector<bool>& used) override {
        auto acc = system.compute_accelerations();

        for (size_t i = 0; i < system.bodies.size(); ++i) {
            if (used[i]) continue;  // Saltar cuerpos marcados
            
            system.bodies[i].velocity = system.bodies[i].velocity + dt * acc[i];
            system.bodies[i].position = system.bodies[i].position + dt * system.bodies[i].velocity;
        }
    }
};

