// regularization/chain/chain_state.h
#pragma once
#include "vec4.h"
#include "vec3.h"

/**
 * @brief Estado de la cadena regularizada para 3 cuerpos.
 * 
 * Contiene:
 * - Variables KS para los dos eslabones (Q1,P1) y (Q2,P2)
 * - Energía interna del sistema de 3 cuerpos
 * - Tiempo ficticio tau
 * - Centro de masa y su tiempo físico
 * - Masas de los tres cuerpos (necesarias para las derivadas)
 */
struct Chain3State {
    // Variables KS para los dos eslabones
    Vec4 Q1, P1;  // Primer eslabón (entre cuerpo 1 y 2)
    Vec4 Q2, P2;  // Segundo eslabón (entre cuerpo 2 y 3)

    // Energía interna del sistema de 3 cuerpos (variable)
    double energy;

    // Tiempo ficticio
    double tau;

    // Centro de masa y su tiempo físico
    Vec3 cm_pos;
    Vec3 cm_vel;
    double cm_time;  // Tiempo físico para el CM (sincronización)

    // Masas de los tres cuerpos (en orden de la cadena)
    double masses[3];

    // Constructor por defecto
    Chain3State() 
        : Q1(), P1(), Q2(), P2()
        , energy(0.0), tau(0.0)
        , cm_pos(), cm_vel(), cm_time(0.0)
    {
        masses[0] = masses[1] = masses[2] = 0.0;
    }

    // Constructor con parámetros básicos
    Chain3State(const Vec4& q1, const Vec4& p1, const Vec4& q2, const Vec4& p2,
                double e, double t, const Vec3& cm_p, const Vec3& cm_v, double cm_t,
                double m1, double m2, double m3)
        : Q1(q1), P1(p1), Q2(q2), P2(p2)
        , energy(e), tau(t)
        , cm_pos(cm_p), cm_vel(cm_v), cm_time(cm_t)
    {
        masses[0] = m1;
        masses[1] = m2;
        masses[2] = m3;
    }

    // Métodos de acceso a masas
    double m1() const { return masses[0]; }
    double m2() const { return masses[1]; }
    double m3() const { return masses[2]; }
    double total_mass() const { return masses[0] + masses[1] + masses[2]; }

    // Validación básica
    bool is_valid() const {
        return (m1() > 0.0 && m2() > 0.0 && m3() > 0.0);
    }
};