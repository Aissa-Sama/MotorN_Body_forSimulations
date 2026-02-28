// integrators/leapfrog_integrator.cpp
#include "leapfrog_integrator.h"
#include "nbody_system.h"

void LeapfrogIntegrator::step(
    NBodySystem& system,
    double dt,
    const std::vector<bool>& used
) {
    const size_t N = system.bodies.size();
    system.invalidate_accelerations();
    auto acc = system.compute_accelerations();
    for (size_t i = 0; i < N; ++i) {
        if (used[i]) continue;
        system.bodies[i].velocity += (dt * 0.5) * acc[i];
    }
    for (size_t i = 0; i < N; ++i) {
        if (used[i]) continue;
        system.bodies[i].position += dt * system.bodies[i].velocity;
    }
    system.invalidate_accelerations();
    acc = system.compute_accelerations();
    for (size_t i = 0; i < N; ++i) {
        if (used[i]) continue;
        system.bodies[i].velocity += (dt * 0.5) * acc[i];
    }
}