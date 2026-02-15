// integrators/rk45_integrator.h
#pragma once
#include <vector>
#include "integrator.h"
#include "nbody_system.h"

class RK45Integrator : public Integrator {
public:
    struct Parameters {
        double abs_tol = 1e-9;
        double rel_tol = 1e-9;
        double safety = 0.9;
        double min_dt = 1e-8;
        double max_dt = 0.1;
        double min_scale = 0.2;
        double max_scale = 5.0;
    };

    explicit RK45Integrator(Parameters params = Parameters{});
    
    void step(NBodySystem& system, double dt, const std::vector<bool>& used) override;

private:
    struct State {
        std::vector<Vec3> pos;
        std::vector<Vec3> vel;
    };
    
    struct Derivative {
        std::vector<Vec3> dpos_dt; // = velocity
        std::vector<Vec3> dvel_dt; // = acceleration
    };
    
    Parameters params;
    
    State capture_state(const NBodySystem& system, const std::vector<bool>& used) const;
    void restore_state(NBodySystem& system, const State& s, const std::vector<bool>& used) const;
    
    Derivative compute_derivatives(const NBodySystem& system, const State& current, 
                                   const std::vector<bool>& used) const;
    
    // Coeficientes de Dormand-Prince (RK5(4)7M)
    static constexpr double a2 = 1.0/5.0;
    static constexpr double a3 = 3.0/40.0;
    static constexpr double a4 = 44.0/45.0;
    static constexpr double a5 = 19372.0/6561.0;
    static constexpr double a6 = 9017.0/3168.0;
    static constexpr double a7 = 35.0/384.0;
    
    static constexpr double b21 = 1.0/5.0;
    
    static constexpr double b31 = 3.0/40.0;
    static constexpr double b32 = 9.0/40.0;
    
    static constexpr double b41 = 44.0/45.0;
    static constexpr double b42 = -56.0/15.0;
    static constexpr double b43 = 32.0/9.0;
    
    static constexpr double b51 = 19372.0/6561.0;
    static constexpr double b52 = -25360.0/2187.0;
    static constexpr double b53 = 64448.0/6561.0;
    static constexpr double b54 = -212.0/729.0;
    
    static constexpr double b61 = 9017.0/3168.0;
    static constexpr double b62 = -355.0/33.0;
    static constexpr double b63 = 46732.0/5247.0;
    static constexpr double b64 = 49.0/176.0;
    static constexpr double b65 = -5103.0/18656.0;
    
    static constexpr double b71 = 35.0/384.0;
    static constexpr double b72 = 0.0;
    static constexpr double b73 = 500.0/1113.0;
    static constexpr double b74 = 125.0/192.0;
    static constexpr double b75 = -2187.0/6784.0;
    static constexpr double b76 = 11.0/84.0;
    
    // Coeficientes para el error (diferencia entre orden 4 y 5)
    static constexpr double e1 = 71.0/57600.0 - 35.0/384.0;
    static constexpr double e3 = -71.0/16695.0 - 500.0/1113.0;
    static constexpr double e4 = 71.0/1920.0 - 125.0/192.0;
    static constexpr double e5 = -17253.0/339200.0 + 2187.0/6784.0;
    static constexpr double e6 = 22.0/525.0 - 11.0/84.0;
    static constexpr double e7 = -1.0/40.0;
};