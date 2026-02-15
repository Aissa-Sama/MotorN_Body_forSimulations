#pragma once
#include <vector>
#include <utility>
#include "nbody_system.h"

class BinaryDetector {
public:
    struct Binary {
        size_t i, j;
        double binding_energy;
    };

    std::vector<Binary> detect(const NBodySystem& system) const;

private:
    double relative_energy(
        const Body& a,
        const Body& b,
        double G
    ) const;
};
