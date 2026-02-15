// initial_conditions.h (actualizado)
#pragma once
#include "nbody_system.h"
#include <vector>

class InitialConditions {
public:
    // Sistema solar simple (Sol + planetas)
    static NBodySystem solar_system();
    
    // Sistema binario + estrella de campo
    static NBodySystem binary_with_field();
    
    // Cúmulo Plummer (distribución realista de estrellas)
    static NBodySystem plummer_cluster(int n_bodies);
    
    // Sistema de prueba: 3 cuerpos en figura-8 (órbita estable conocida)
    static NBodySystem figure_eight();
    
    // Sistema aleatorio para testing
    static NBodySystem random_system(int n_bodies, double radius);
    
    // Binaria Kepleriana pura (NUEVO)
    static NBodySystem kepler_binary(double a = 1.0, double e = 0.0);
};