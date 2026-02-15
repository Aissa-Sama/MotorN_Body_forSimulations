// regularization/ks/ks_state.h
#pragma once
#include "vec3.h"

struct KSState {
    // coordenadas KS (4D)
    double u[4];
    double w[4];   // velocidades conjugadas
    
    double energy; // energía relativa constante
    double tau;    // tiempo ficticio
};
