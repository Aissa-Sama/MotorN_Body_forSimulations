#pragma once

#include "integrator.h"

class VelocityVerletIntegrator : public Integrator {
public:
    void step(
        NBodySystem& system,
        double dt,
        const std::vector<bool>& used
    ) override;
};
