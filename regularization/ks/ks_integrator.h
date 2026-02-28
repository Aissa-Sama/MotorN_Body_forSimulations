// regularization/ks/ks_integrator.h
#pragma once
#include "binary_state.h"
#include "ks_state.h"

class KSIntegrator {
public:
    explicit KSIntegrator(double internal_dt);
    void integrate(BinaryState& binary, double dt_global);

protected:
    double internal_dt;
    
    KSState to_ks(const BinaryState& b) const;
    void from_ks(const KSState& ks, BinaryState& b) const;
};