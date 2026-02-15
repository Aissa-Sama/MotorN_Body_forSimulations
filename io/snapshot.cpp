// io/snapshot.cpp
#include "snapshot.h"
#include <fstream>
#include <iomanip>

void SnapshotManager::save(const NBodySystem& system, 
                           const std::string& filename,
                           double time) {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    
    // Formato simple: tiempo, número de cuerpos
    file << time << "\n";
    file << system.bodies.size() << "\n";
    
    // Cada cuerpo: masa, pos.x, pos.y, pos.z, vel.x, vel.y, vel.z
    for (const auto& b : system.bodies) {
        file << std::setprecision(15)
             << b.mass << " "
             << b.position.x << " " << b.position.y << " " << b.position.z << " "
             << b.velocity.x << " " << b.velocity.y << " " << b.velocity.z << "\n";
    }
}

NBodySystem SnapshotManager::load(const std::string& filename, double& time) {
    NBodySystem system;
    std::ifstream file(filename);
    if (!file.is_open()) return system;
    
    file >> time;
    size_t n_bodies;
    file >> n_bodies;
    
    system.bodies.resize(n_bodies);
    for (auto& b : system.bodies) {
        file >> b.mass
             >> b.position.x >> b.position.y >> b.position.z
             >> b.velocity.x >> b.velocity.y >> b.velocity.z;
    }
    
    return system;
}

void SnapshotManager::set_autosave(int interval, const std::string& basename) {
    // Implementación simple: guardaríamos el intervalo y nombre base
    // Por ahora solo almacenamos, la lógica se llamaría desde main
    (void)interval;  // Evitar warning
    (void)basename;
}