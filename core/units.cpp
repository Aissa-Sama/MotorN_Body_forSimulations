// core/units.cpp
#include "units.h"
#include "nbody_system.h"  // <--- ESTE INCLUDE FALTABA
#include "body.h"          // <--- ESTE INCLUDE FALTABA

namespace units {
    NBodySystem to_simulation_units(const NBodySystem& real_system, const NBodyUnits& units) {
        NBodySystem sim_system;
        sim_system.bodies.reserve(real_system.bodies.size());
        
        for (const auto& b : real_system.bodies) {
            Body sim_body;
            sim_body.mass = b.mass / units.mass_unit;
            sim_body.position = b.position / units.length_unit;
            sim_body.velocity = b.velocity / (units.length_unit / units.time_unit);
            sim_system.bodies.push_back(sim_body);
        }
        
        sim_system.G = 1.0;
        
        return sim_system;
    }
}