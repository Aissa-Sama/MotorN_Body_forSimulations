#pragma once

#include <vector>

class NBodySystem;

class Integrator {
public:
    virtual ~Integrator() = default;

    // Máscara: used[i] == true → NO integrar ese cuerpo
    virtual void step(
        NBodySystem& system,
        double dt,
        const std::vector<bool>& used
    ) = 0;
};
