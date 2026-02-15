// regularization/ks/ks_integrator.h (actualizado)
#pragma once
#include "binary_state.h"
#include "ks_state.h"

class KSIntegrator {
public:
    explicit KSIntegrator(double internal_dt);
    
    // Integra binary desde t hasta t+dt_global
    // Usa internal_dt internamente
    void integrate(BinaryState& binary, double dt_global);

private:
    double internal_dt;
    
    KSState to_ks(const BinaryState& b) const;
    void from_ks(const KSState& ks, BinaryState& b) const;
    void step_ks(KSState& ks);  // un paso con internal_dt
    
    // Relación tiempo físico ↔ ficticio
    double physical_time_step(double r, double dt_global) const;
};