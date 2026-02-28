#pragma once
#include <vector>
#include "integrator.h"
#include "nbody_system.h"

class LeapfrogIntegrator : public Integrator {
public:
    void step(
        NBodySystem& system,
        double dt,
        const std::vector<bool>& used
    ) override;
};