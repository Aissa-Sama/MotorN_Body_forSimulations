// integrators/velocity_verlet_integrator.cpp
#include "velocity_verlet_integrator.h"
#include "nbody_system.h"

void VelocityVerletIntegrator::step(
    NBodySystem& system,
    double dt,
    const std::vector<bool>& used
) {
    const int N = static_cast<int>(system.bodies.size());
    
    // Obtener aceleraciones actuales
    auto acc = system.compute_accelerations();
    
    // --- Actualizar posiciones (usando aceleraciones actuales) ---
    for (int i = 0; i < N; ++i) {
        if (used[i]) continue;
        
        auto& b = system.bodies[i];
        b.position = b.position + b.velocity * dt + 0.5 * acc[i] * dt * dt;
    }
    
    // Invalidar caché porque las posiciones cambiaron
    system.invalidate_accelerations();
    
    // Calcular nuevas aceleraciones
    auto new_acc = system.compute_accelerations();
    
    // --- Actualizar velocidades (usando nuevas aceleraciones) ---
    for (int i = 0; i < N; ++i) {
        if (used[i]) continue;
        
        auto& b = system.bodies[i];
        b.velocity = b.velocity + 0.5 * (acc[i] + new_acc[i]) * dt;
    }
}
