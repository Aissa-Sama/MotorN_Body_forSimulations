// integrators/rk4_integrator.h
#pragma once
#include "integrator.h"

class RK4Integrator : public Integrator {
public:
    void step(NBodySystem& system, double dt, const std::vector<bool>& used) override {
        const size_t N = system.bodies.size();

        // Estado original (solo de cuerpos no usados)
        std::vector<Vec3> x0(N), v0(N);
        std::vector<size_t> active_indices;
        for (size_t i = 0; i < N; ++i) {
            if (!used[i]) {
                x0[i] = system.bodies[i].position;
                v0[i] = system.bodies[i].velocity;
                active_indices.push_back(i);
            }
        }

        // k1 (usando aceleraciones actuales)
        auto a1 = system.compute_accelerations();
        std::vector<Vec3> k1x(N), k1v(N);
        for (size_t i = 0; i < N; ++i) {
            if (!used[i]) {
                k1x[i] = v0[i];
                k1v[i] = a1[i];
            }
        }

        // k2
        apply_temp_state(system, x0, v0, k1x, k1v, dt * 0.5, used);
        auto a2 = system.compute_accelerations();
        std::vector<Vec3> k2x(N), k2v(N);
        for (size_t i = 0; i < N; ++i) {
            if (!used[i]) {
                k2x[i] = system.bodies[i].velocity;
                k2v[i] = a2[i];
            }
        }

        // k3
        apply_temp_state(system, x0, v0, k2x, k2v, dt * 0.5, used);
        auto a3 = system.compute_accelerations();
        std::vector<Vec3> k3x(N), k3v(N);
        for (size_t i = 0; i < N; ++i) {
            if (!used[i]) {
                k3x[i] = system.bodies[i].velocity;
                k3v[i] = a3[i];
            }
        }

        // k4
        apply_temp_state(system, x0, v0, k3x, k3v, dt, used);
        auto a4 = system.compute_accelerations();
        std::vector<Vec3> k4x(N), k4v(N);
        for (size_t i = 0; i < N; ++i) {
            if (!used[i]) {
                k4x[i] = system.bodies[i].velocity;
                k4v[i] = a4[i];
            }
        }

        // Estado final (solo para cuerpos activos)
        for (size_t i = 0; i < N; ++i) {
            if (!used[i]) {
                system.bodies[i].position = x0[i] + (dt / 6.0) * 
                    (k1x[i] + 2.0*k2x[i] + 2.0*k3x[i] + k4x[i]);

                system.bodies[i].velocity = v0[i] + (dt / 6.0) * 
                    (k1v[i] + 2.0*k2v[i] + 2.0*k3v[i] + k4v[i]);
            }
        }
    }

private:
    void apply_temp_state(
        NBodySystem& system,
        const std::vector<Vec3>& x0,
        const std::vector<Vec3>& v0,
        const std::vector<Vec3>& kx,
        const std::vector<Vec3>& kv,
        double dt,
        const std::vector<bool>& used)
    {
        for (size_t i = 0; i < system.bodies.size(); ++i) {
            if (!used[i]) {
                system.bodies[i].position = x0[i] + dt * kx[i];
                system.bodies[i].velocity = v0[i] + dt * kv[i];
            }
        }
    }
};