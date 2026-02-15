// integrators/rk45_integrator.cpp
#include "rk45_integrator.h"
#include <cmath>
#include <algorithm>

RK45Integrator::RK45Integrator(Parameters params)
    : params(params) {}

RK45Integrator::State RK45Integrator::capture_state(
    const NBodySystem& system, 
    const std::vector<bool>& used) const 
{
    State s;
    for (size_t i = 0; i < system.bodies.size(); ++i) {
        if (!used[i]) {
            s.pos.push_back(system.bodies[i].position);
            s.vel.push_back(system.bodies[i].velocity);
        }
    }
    return s;
}

void RK45Integrator::restore_state(
    NBodySystem& system, 
    const State& s, 
    const std::vector<bool>& used) const 
{
    size_t idx = 0;
    for (size_t i = 0; i < system.bodies.size(); ++i) {
        if (!used[i]) {
            system.bodies[i].position = s.pos[idx];
            system.bodies[i].velocity = s.vel[idx];
            idx++;
        }
    }
}

RK45Integrator::Derivative RK45Integrator::compute_derivatives(
    const NBodySystem& system, 
    const State& current,
    const std::vector<bool>& used) const 
{
    Derivative d;
    size_t n_active = current.pos.size();
    d.dpos_dt.reserve(n_active);
    d.dvel_dt.reserve(n_active);
    
    // Las velocidades son las derivadas de la posición
    for (size_t i = 0; i < n_active; ++i) {
        d.dpos_dt.push_back(current.vel[i]);
    }
    
    // Calcular aceleraciones
    // Necesitamos reconstruir el sistema completo para las aceleraciones
    NBodySystem temp;
    temp.G = system.G;
    
    // Primero, copiar todos los cuerpos
    temp.bodies = system.bodies;
    
    // Luego, actualizar los activos con current
    size_t idx = 0;
    for (size_t i = 0; i < temp.bodies.size(); ++i) {
        if (!used[i]) {
            temp.bodies[i].position = current.pos[idx];
            temp.bodies[i].velocity = current.vel[idx];
            idx++;
        }
    }
    
    // Calcular aceleraciones
    auto acc = temp.compute_accelerations();
    
    // Extraer solo las aceleraciones de los activos
    idx = 0;
    for (size_t i = 0; i < temp.bodies.size(); ++i) {
        if (!used[i]) {
            d.dvel_dt.push_back(acc[i]);
            idx++;
        }
    }
    
    return d;
}

void RK45Integrator::step(NBodySystem& system, double dt, const std::vector<bool>& used) {
    // Si dt es 0, no hacemos nada
    if (dt == 0.0) return;
    
    State y0 = capture_state(system, used);
    size_t n_active = y0.pos.size();
    if (n_active == 0) return;
    
    double current_dt = dt;
    double time_done = 0.0;
    
    while (time_done < dt - 1e-12) {
        // Asegurar que no nos pasamos
        current_dt = std::min(current_dt, dt - time_done);
        
        // Etapa 1: k1
        Derivative k1 = compute_derivatives(system, y0, used);
        
        // Etapa 2
        State y2;
        for (size_t i = 0; i < n_active; ++i) {
            Vec3 pos = y0.pos[i] + current_dt * b21 * k1.dpos_dt[i];
            Vec3 vel = y0.vel[i] + current_dt * b21 * k1.dvel_dt[i];
            y2.pos.push_back(pos);
            y2.vel.push_back(vel);
        }
        Derivative k2 = compute_derivatives(system, y2, used);
        
        // Etapa 3
        State y3;
        for (size_t i = 0; i < n_active; ++i) {
            Vec3 pos = y0.pos[i] + current_dt * (b31 * k1.dpos_dt[i] + b32 * k2.dpos_dt[i]);
            Vec3 vel = y0.vel[i] + current_dt * (b31 * k1.dvel_dt[i] + b32 * k2.dvel_dt[i]);
            y3.pos.push_back(pos);
            y3.vel.push_back(vel);
        }
        Derivative k3 = compute_derivatives(system, y3, used);
        
        // Etapa 4
        State y4;
        for (size_t i = 0; i < n_active; ++i) {
            Vec3 pos = y0.pos[i] + current_dt * (b41 * k1.dpos_dt[i] + b42 * k2.dpos_dt[i] + 
                                                  b43 * k3.dpos_dt[i]);
            Vec3 vel = y0.vel[i] + current_dt * (b41 * k1.dvel_dt[i] + b42 * k2.dvel_dt[i] + 
                                                  b43 * k3.dvel_dt[i]);
            y4.pos.push_back(pos);
            y4.vel.push_back(vel);
        }
        Derivative k4 = compute_derivatives(system, y4, used);
        
        // Etapa 5
        State y5;
        for (size_t i = 0; i < n_active; ++i) {
            Vec3 pos = y0.pos[i] + current_dt * (b51 * k1.dpos_dt[i] + b52 * k2.dpos_dt[i] + 
                                                  b53 * k3.dpos_dt[i] + b54 * k4.dpos_dt[i]);
            Vec3 vel = y0.vel[i] + current_dt * (b51 * k1.dvel_dt[i] + b52 * k2.dvel_dt[i] + 
                                                  b53 * k3.dvel_dt[i] + b54 * k4.dvel_dt[i]);
            y5.pos.push_back(pos);
            y5.vel.push_back(vel);
        }
        Derivative k5 = compute_derivatives(system, y5, used);
        
        // Etapa 6
        State y6;
        for (size_t i = 0; i < n_active; ++i) {
            Vec3 pos = y0.pos[i] + current_dt * (b61 * k1.dpos_dt[i] + b62 * k2.dpos_dt[i] + 
                                                  b63 * k3.dpos_dt[i] + b64 * k4.dpos_dt[i] + 
                                                  b65 * k5.dpos_dt[i]);
            Vec3 vel = y0.vel[i] + current_dt * (b61 * k1.dvel_dt[i] + b62 * k2.dvel_dt[i] + 
                                                  b63 * k3.dvel_dt[i] + b64 * k4.dvel_dt[i] + 
                                                  b65 * k5.dvel_dt[i]);
            y6.pos.push_back(pos);
            y6.vel.push_back(vel);
        }
        Derivative k6 = compute_derivatives(system, y6, used);
        
        // Etapa 7 (solo para el cálculo del error)
        State y7;
        for (size_t i = 0; i < n_active; ++i) {
            Vec3 pos = y0.pos[i] + current_dt * (b71 * k1.dpos_dt[i] + b73 * k3.dpos_dt[i] + 
                                                  b74 * k4.dpos_dt[i] + b75 * k5.dpos_dt[i] + 
                                                  b76 * k6.dpos_dt[i]);
            Vec3 vel = y0.vel[i] + current_dt * (b71 * k1.dvel_dt[i] + b73 * k3.dvel_dt[i] + 
                                                  b74 * k4.dvel_dt[i] + b75 * k5.dvel_dt[i] + 
                                                  b76 * k6.dvel_dt[i]);
            y7.pos.push_back(pos);
            y7.vel.push_back(vel);
        }
        Derivative k7 = compute_derivatives(system, y7, used);
        // Calcular error (diferencia entre orden 5 y orden 4)
        double max_error = 0.0;
        for (size_t i = 0; i < n_active; ++i) {
            // Estimación del error usando los coeficientes e
            Vec3 error_pos = (e1 * k1.dpos_dt[i] + e3 * k3.dpos_dt[i] + 
                              e4 * k4.dpos_dt[i] + e5 * k5.dpos_dt[i] + 
                              e6 * k6.dpos_dt[i] + e7 * k7.dpos_dt[i]) * current_dt;
            
            double scale = params.abs_tol + params.rel_tol * 
                          std::max(norm(y0.pos[i]), norm(y7.pos[i]));
            double err = norm(error_pos) / scale;
            max_error = std::max(max_error, err);
        }
        
        // Aceptar o rechazar el paso
        if (max_error <= 1.0 || current_dt <= params.min_dt) {
            // Paso aceptado
            y0 = y7;
            time_done += current_dt;
        }
        
        // Calcular nuevo dt
        if (max_error > 0) {
            double scale = params.safety * std::pow(1.0 / max_error, 0.2);
            scale = std::clamp(scale, params.min_scale, params.max_scale);
            current_dt *= scale;
            current_dt = std::clamp(current_dt, params.min_dt, params.max_dt);
        }
    }
    
    // Escribir resultado final
    restore_state(system, y0, used);
}
    