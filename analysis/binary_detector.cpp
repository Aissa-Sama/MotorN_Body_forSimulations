#include "binary_detector.h"
#include "vec3.h"

double BinaryDetector::relative_energy(
    const Body& a,
    const Body& b,
    double G
) const {
    Vec3 r = b.position - a.position;
    Vec3 v = b.velocity - a.velocity;

    double mu = (a.mass * b.mass) / (a.mass + b.mass);
    double kinetic = 0.5 * mu * dot(v, v);
    double potential = -G * a.mass * b.mass / norm(r);

    return kinetic + potential;
}

std::vector<BinaryDetector::Binary>
BinaryDetector::detect(const NBodySystem& system) const {
    std::vector<Binary> binaries;

    for (size_t i = 0; i < system.bodies.size(); ++i) {
        for (size_t j = i + 1; j < system.bodies.size(); ++j) {
            double E = relative_energy(
                system.bodies[i],
                system.bodies[j],
                system.G
            );

            if (E < 0.0) {
                binaries.push_back({i, j, E});
            }
        }
    }
    return binaries;
}
