// core/units.h
#pragma once
#include "nbody_system.h"  // <--- ESTE INCLUDE FALTABA
#include <string>

namespace units {
    // Constantes físicas
    constexpr double G_SI = 6.67430e-11;  // m^3 kg^-1 s^-2
    
    // Unidades astronómicas
    constexpr double AU = 1.495978707e11;  // metros
    constexpr double solar_mass = 1.9885e30;  // kg
    constexpr double year = 31557600.0;  // segundos
    
    // Sistema de unidades adaptado (G=1, masas=1, distancias=1)
    struct NBodyUnits {
        double mass_unit;    // ej: masa solar
        double length_unit;  // ej: AU
        double time_unit;    // ej: año
        double G = 1.0;      // en estas unidades
    };
    
    // Convertir sistema a unidades de simulación
    NBodySystem to_simulation_units(const NBodySystem& real_system, 
                                     const NBodyUnits& units);
}