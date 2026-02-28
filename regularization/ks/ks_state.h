// regularization/ks/ks_state.h
#pragma once
#include "vec3.h"

struct KSState {
    // coordenadas KS (4D)
    double u[4];
    double w[4];   // velocidades conjugadas KS
    
    // CONVENCIÓN: energy se almacena NEGADA respecto a la energía física.
    //   energy_stored = -E_física > 0  para órbitas ligadas (E_física < 0)
    // Esto hace que el leapfrog del oscilador armónico sea:
    //   w[j] -= ½·dτ·energy·u[j]
    // con energy > 0, lo que produce confinamiento (oscilador estable).
    // Al añadir perturbaciones, tener en cuenta este signo.
    double energy;

    double tau;    // tiempo ficticio KS
};

