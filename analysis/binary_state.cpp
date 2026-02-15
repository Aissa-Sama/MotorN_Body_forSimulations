// analysis/binary_state.cpp
#include "binary_state.h"

BinaryState::BinaryState(const Body& a, const Body& b) {
    m1 = a.mass;
    m2 = b.mass;

    r = b.position - a.position;
    v = b.velocity - a.velocity;

    double M = m1 + m2;
    R = (m1 * a.position + m2 * b.position) / M;
    V = (m1 * a.velocity + m2 * b.velocity) / M;
}

void BinaryState::integrate(double /*dt*/) {
    // Placeholder: la dinámica real la hace KSIntegrator
}

void BinaryState::write_back(Body& a, Body& b) const {
    double M = m1 + m2;

    a.position = R - (m2 / M) * r;
    b.position = R + (m1 / M) * r;

    a.velocity = V - (m2 / M) * v;
    b.velocity = V + (m1 / M) * v;
}

double BinaryState::separation() const {
    return norm(r);
}

// Nota: advance_cm ahora está implementado inline en el header
