// analysis/binary_state.h
#pragma once
#include "body.h"
#include "vec3.h"

class BinaryState {
public:
    BinaryState(const Body& a, const Body& b);
    
    // Integración (usado por KSIntegrator)
    void integrate(double dt);
    void write_back(Body& a, Body& b) const;
    
    // Getters
    double separation() const;
    
    // Getters para masas
    double mass1() const { return m1; }
    double mass2() const { return m2; }
    double total_mass() const { return m1 + m2; }
    double reduced_mass() const { return (m1 * m2) / (m1 + m2); }
    
    // Getters para posición y velocidad relativas
    const Vec3& relative_position() const { return r; }
    const Vec3& relative_velocity() const { return v; }
    
    // Getters para centro de masa
    const Vec3& center_of_mass() const { return R; }
    const Vec3& cm_velocity() const { return V; }
    
    // Setters para KS (necesarios para modificar r y v)
    void set_relative_position(const Vec3& new_r) { r = new_r; }
    void set_relative_velocity(const Vec3& new_v) { v = new_v; }
    void set_center_of_mass(const Vec3& new_R) { R = new_R; }
    void set_cm_velocity(const Vec3& new_V) { V = new_V; }
    
    // Actualizar centro de masa (movimiento uniforme)
    void advance_cm(double dt) {
        R = R + V * dt;
    }

private:
    double m1, m2;
    Vec3 r;   // posición relativa
    Vec3 v;   // velocidad relativa
    Vec3 R;   // centro de masa
    Vec3 V;   // velocidad del CM
};
