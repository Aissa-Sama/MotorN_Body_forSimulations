// analysis/system_analyzer.h - VERSIÓN CORREGIDA
#pragma once
#include "nbody_system.h"
#include "binary_detector.h"
#include <vector>

class SystemAnalyzer {
public:
    struct Snapshot {
        double time;
        double energy;
        Vec3 momentum;           // <-- Necesita #include "vec3.h"
        Vec3 angular_momentum;    // <-- Necesita #include "vec3.h"
        double virial_ratio;
        int n_binaries;
        double min_separation;
        double max_separation;
    };
    
    // Analizar sistema completo
    Snapshot analyze(const NBodySystem& system, double time);
    
    // Detectar binarias (usando el detector que ya tienes)
    std::vector<BinaryDetector::Binary> find_binaries(const NBodySystem& system);
    
    // Calcular radio de medio masa
    double half_mass_radius(const NBodySystem& system);
};