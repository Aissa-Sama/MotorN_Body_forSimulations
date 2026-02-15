// analysis/system_analyzer.cpp - VERSIÓN CORREGIDA
#include "system_analyzer.h"
#include "binary_detector.h"
#include "vec3.h"           
#include <algorithm>        
#include <vector>
#include <cmath>

SystemAnalyzer::Snapshot SystemAnalyzer::analyze(const NBodySystem& system, double time) {
    Snapshot s;
    s.time = time;
    s.energy = system.total_energy();
    s.momentum = system.total_momentum();
    s.angular_momentum = system.total_angular_momentum();
    
    double T = system.kinetic_energy();
    double V = system.potential_energy();
    s.virial_ratio = (V != 0.0) ? -2.0 * T / V : 0.0;
    
    auto binaries = find_binaries(system);
    s.n_binaries = static_cast<int>(binaries.size());
    
    s.min_separation = 1e99;
    s.max_separation = 0.0;
    for (const auto& b : binaries) {
        Vec3 r = system.bodies[b.j].position - system.bodies[b.i].position;
        double d = norm(r);
        s.min_separation = std::min(s.min_separation, d);
        s.max_separation = std::max(s.max_separation, d);
    }
    
    return s;
}

// IMPLEMENTACIÓN QUE FALTABA
std::vector<BinaryDetector::Binary> SystemAnalyzer::find_binaries(const NBodySystem& system) {
    BinaryDetector detector;
    return detector.detect(system);
}

double SystemAnalyzer::half_mass_radius(const NBodySystem& system) {
    if (system.bodies.empty()) return 0.0;
    
    double total_mass = 0.0;
    for (const auto& b : system.bodies) {
        total_mass += b.mass;
    }
    
    // Calcular centro de masa
    Vec3 cm = {0, 0, 0};
    for (const auto& b : system.bodies) {
        cm = cm + b.mass * b.position;
    }
    cm = cm / total_mass;
    
    std::vector<double> distances;
    for (const auto& b : system.bodies) {
        distances.push_back(norm(b.position - cm));
    }
    std::sort(distances.begin(), distances.end());
    
    double half_mass = total_mass * 0.5;
    double cumulative = 0.0;
    double mass_per_body = total_mass / system.bodies.size();
    
    for (double d : distances) {
        cumulative += mass_per_body;
        if (cumulative >= half_mass) return d;
    }
    
    return distances.back();
}