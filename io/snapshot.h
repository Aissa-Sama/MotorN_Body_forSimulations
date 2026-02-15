// io/snapshot.h
#pragma once
#include "nbody_system.h"  // <--- ESTE INCLUDE FALTABA
#include <string>          // <--- ESTE INCLUDE FALTABA

class SnapshotManager {
public:
    // Guardar estado actual
    void save(const NBodySystem& system, 
              const std::string& filename,
              double time);
    
    // Cargar estado previo
    NBodySystem load(const std::string& filename, double& time);
    
    // Guardar cada N pasos (autosave)
    void set_autosave(int interval, const std::string& basename);
};