// core/nbody_system.cpp - IMPLEMENTACIÓN COMPLETA
#include "nbody_system.h"
#include <cmath>

std::vector<Vec3> NBodySystem::compute_accelerations() const {
    const size_t N = bodies.size();
    if (accelerations_valid && accelerations.size() == N) {
        return accelerations;
    }
    accelerations.assign(N, {0, 0, 0});
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (i == j) continue;
            Vec3 r = bodies[j].position - bodies[i].position;
            double d = r.norm();
            double d3 = d * d * d;
            accelerations[i] = accelerations[i] + (G * bodies[j].mass / d3) * r;
        }
    }
    accelerations_valid = true;
    return accelerations;
}

void NBodySystem::invalidate_accelerations() {
    accelerations_valid = false;
}

double NBodySystem::kinetic_energy() const {
    double T = 0.0;
    for (const auto& b : bodies) {
        T += 0.5 * b.mass * b.velocity.norm2();
    }
    return T;
}

double NBodySystem::potential_energy() const {
    double V = 0.0;
    const size_t N = bodies.size();
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = i + 1; j < N; ++j) {
            Vec3 r = bodies[j].position - bodies[i].position;
            double d = r.norm();
            V -= G * bodies[i].mass * bodies[j].mass / d;
        }
    }
    return V;
}

double NBodySystem::total_energy() const {
    return kinetic_energy() + potential_energy();
}

Vec3 NBodySystem::total_momentum() const {
    Vec3 P = {0, 0, 0};
    for (const auto& b : bodies) {
        P = P + b.mass * b.velocity;
    }
    return P;
}

Vec3 NBodySystem::total_angular_momentum() const {
    Vec3 L = {0, 0, 0};
    for (const auto& b : bodies) {
        Vec3 r = b.position;
        Vec3 v = b.velocity;
        Vec3 rxv = {
            r.y * v.z - r.z * v.y,
            r.z * v.x - r.x * v.z,
            r.x * v.y - r.y * v.x
        };
        L = L + b.mass * rxv;
    }
    return L;
}