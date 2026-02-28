#pragma once

#include <vector>
#include "nbody_system.h"

class Integrator {
public:
    virtual void step(
        NBodySystem& system,
        double dt,
        const std::vector<bool>& used
    ) = 0;

    virtual ~Integrator() = default;
};