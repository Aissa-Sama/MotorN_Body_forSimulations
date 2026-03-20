// tests/test_archain_n_ks.cpp
// Fase 6B — KS dentro del AR-chain. Criterios calibrados con análisis del splitting.
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>
#include <chrono>
#include "archain_n_ks_integrator.h"
#include "archain_n_bs_integrator.h"
#include "nbody_system.h"
#include "body.h"

static NBodySystem make_binary_ecc(double m1, double m2, double a, double e) {
    const double M = m1 + m2;
    const double r_apo = a*(1+e);
    const double v_apo = std::sqrt(M*(1-e)/(a*(1+e)));
    NBodySystem sys; sys.bodies.resize(2);
    sys.bodies[0] = {Vec3( m2/M*r_apo, 0,0), Vec3(0, m2/M*v_apo, 0), m1};
    sys.bodies[1] = {Vec3(-m1/M*r_apo, 0,0), Vec3(0,-m1/M*v_apo, 0), m2};
    return sys;
}

static NBodySystem make_triple_compact(double m1,double m2,double m3,
                                        double a_in,double e_in,double a_out) {
    const double M12=m1+m2, Mtot=M12+m3;
    const double r_apo=a_in*(1+e_in);
    const double v_apo=std::sqrt(M12*(1-e_in)/(a_in*(1+e_in)));
    const double v3=std::sqrt(Mtot/a_out)*0.6;
    NBodySystem sys; sys.bodies.resize(3);
    sys.bodies[0]={Vec3( m2/M12*r_apo,0,0),Vec3(0, m2/M12*v_apo,0),m1};
    sys.bodies[1]={Vec3(-m1/M12*r_apo,0,0),Vec3(0,-m1/M12*v_apo,0),m2};
    sys.bodies[2]={Vec3(a_out,0,0),Vec3(0,v3,0),m3};
    return sys;
}

static double compute_energy(const NBodySystem& sys) {
    double T=0,U=0; const int N=sys.bodies.size();
    for(int i=0;i<N;++i){ T+=0.5*sys.bodies[i].mass*dot(sys.bodies[i].velocity,sys.bodies[i].velocity);
        for(int j=i+1;j<N;++j){Vec3 dr=sys.bodies[j].position-sys.bodies[i].position;U-=sys.G*sys.bodies[i].mass*sys.bodies[j].mass/norm(dr);}}
    return T+U;
}
static bool PASS(const char* n){std::cout<<"  [PASS] "<<n<<"\n";return true;}
static bool FAIL(const char* n,const std::string& r){std::cout<<"  [FAIL] "<<n<<" — "<<r<<"\n";return false;}

// T1 — Circular (e=0): KS inactivo, energía < 1e-10
static bool test1() {
    const double a=1.0,e=0.0,T=2*M_PI*std::sqrt(a*a*a);
    ARChainNKSIntegrator::KSParameters ks; ks.alpha_ks=0.05;
    ARChainNBSIntegrator::BSParameters bs; bs.bs_eps=1e-12;
    ARChainNKSIntegrator integ(1e-3,bs,ks);
    auto sys=make_binary_ecc(0.5,0.5,a,e);
    const std::vector<int> idx={0,1}; const double E0=compute_energy(sys);
    ARChainNState st=integ.initialize(sys,idx);
    const int n=50; const double dt=T*5/n;
    for(int s=0;s<n;++s) integ.integrate_to_ks(st,st.t_phys+dt);
    integ.write_back(st,sys,idx);
    const double dE=std::abs(compute_energy(sys)-E0)/std::abs(E0);
    std::cout<<"    e=0  |dE/E0|="<<dE<<"  ks_active="<<(integ.ks_was_active()?"true":"false")<<"\n";
    if(integ.ks_was_active()) return FAIL("circular","KS activo para e=0");
    if(dE>1e-10) return FAIL("circular","|dE|>1e-10 en 5 orbitas");
    return PASS("circular (e=0): KS inactivo, |dE/E0|<1e-10");
}

// T2 — Excéntrica (e=0.9): KS se activa, dE<5e-2 con pasos moderados
static bool test2() {
    const double a=1.0,e=0.9,T=2*M_PI*std::sqrt(a*a*a);
    ARChainNKSIntegrator::KSParameters ks; ks.alpha_ks=0.15; ks.hysteresis_factor=2.0; ks.ks_substeps=128;
    ARChainNBSIntegrator::BSParameters bs; bs.bs_eps=1e-11; bs.max_steps=2000000;
    ARChainNKSIntegrator integ(1e-3,bs,ks);
    auto sys=make_binary_ecc(0.5,0.5,a,e);
    const std::vector<int> idx={0,1}; const double E0=compute_energy(sys);
    ARChainNState st=integ.initialize(sys,idx);
    const int n=200; const double dt=T/n; bool activated=false;
    for(int s=0;s<n;++s){integ.integrate_to_ks(st,st.t_phys+dt);if(integ.ks_was_active())activated=true;}
    integ.write_back(st,sys,idx);
    const double dE=std::abs(compute_energy(sys)-E0)/std::abs(E0);
    std::cout<<"    e=0.9  KS_act="<<(activated?"SI":"NO")<<"  |dE/E0|="<<dE
             <<"  (splitting O(dt^2)=O("<<(dt*dt)<<"))\n";
    if(!activated) return FAIL("excentrica","KS no se activo (alpha=0.15 > r_peri/a=0.1)");
    if(dE>5e-2) return FAIL("excentrica","|dE|>5e-2 — splitting demasiado grande");
    return PASS("excentrica (e=0.9): KS activo en perihelio, splitting aceptable");
}

// T3 — Alta excentricidad (e=0.999): dE<1e-3 en 1 período
static bool test3() {
    const double a=1.0,e=0.999,T=2*M_PI*std::sqrt(a*a*a);
    ARChainNKSIntegrator::KSParameters ks; ks.alpha_ks=0.5; ks.hysteresis_factor=2.0; ks.ks_substeps=256;
    ARChainNBSIntegrator::BSParameters bs; bs.bs_eps=1e-10; bs.max_steps=5000000;
    ARChainNKSIntegrator integ(1e-4,bs,ks);
    auto sys=make_binary_ecc(0.5,0.5,a,e);
    const std::vector<int> idx={0,1}; const double E0=compute_energy(sys);
    ARChainNState st=integ.initialize(sys,idx);
    const int n=20; const double dt=T/n;
    for(int s=0;s<n;++s) integ.integrate_to_ks(st,st.t_phys+dt);
    integ.write_back(st,sys,idx);
    const double dE=std::abs(compute_energy(sys)-E0)/std::abs(E0);
    std::cout<<"    e=0.999  r_peri="<<(a*(1-e))<<"  |dE/E0|="<<dE
             <<"  ks="<<(integ.ks_was_active()?"ON":"OFF")<<"\n";
    if(dE>1e-3) return FAIL("e=0.999","|dE|>1e-3 en 1 periodo con KS");
    return PASS("e=0.999: KS activo, |dE/E0|<1e-3 en 1 periodo");
}

// T4 — Triple compacto (e=0.5): KS detecta el par correcto, sin divergencia
//
// e_in=0.5: r_peri=0.025 → ~28K pasos. Caso base que ya funciona.
static bool test4() {
    const double a_in=0.05,e_in=0.5,a_out=1.0;
    const double T_in=2*M_PI*std::sqrt(a_in*a_in*a_in);
    ARChainNKSIntegrator::KSParameters ks; ks.alpha_ks=0.6; ks.hysteresis_factor=2.0; ks.ks_substeps=64;
    ARChainNBSIntegrator::BSParameters bs; bs.bs_eps=1e-10; bs.max_steps=100000;
    ARChainNKSIntegrator integ(1e-3,bs,ks);
    auto sys=make_triple_compact(0.5,0.5,0.2,a_in,e_in,a_out);
    const std::vector<int> idx={0,1,2}; const double E0=compute_energy(sys);
    ARChainNState st=integ.initialize(sys,idx);
    const int n=20; const double dt=T_in/n; bool activated=false;
    for(int s=0;s<n;++s){integ.integrate_to_ks(st,st.t_phys+dt);if(integ.ks_was_active())activated=true;}
    integ.write_back(st,sys,idx);
    const double dE=std::abs(compute_energy(sys)-E0)/std::abs(E0);
    std::cout<<"    a_in="<<a_in<<"  e_in="<<e_in<<"  a_out="<<a_out
             <<"  KS_act="<<(activated?"SI":"NO")
             <<"  link_crit="<<integ.ks_critical_link()
             <<"  |dE/E0|="<<dE<<"\n";
    if(!activated) return FAIL("triple e=0.5","KS no se activo en el par interno");
    if(dE>1e-1) return FAIL("triple e=0.5","|dE|>1e-1 con paso adaptativo KS");
    return PASS("triple (e=0.5): KS activo en par interno, sin divergencia");
}

// T4b — Triple compacto (e=0.9): CASO PROBLEMATICO documentado en Maestro v11
//
// Con e_in=0.9, r_peri=a_in*(1-e)=0.005. Sin el paso adaptativo KS en
// integrate_until(), el AR-chain necesita ~700K pasos por periodo (~15 min).
// Con el paso adaptativo KS activo en el bucle GBS (integrate_to_ks), el
// numero de pasos cae a O(1/sqrt(r_peri)) — objetivo: < 60s.
//
// Criterios:
//   - KS se activa para el par interno         (robustez)
//   - |dE/E0| < 1e-1                           (no divergencia)
//   - Tiempo de ejecucion < 60s                (rendimiento — la deuda tecnica)
static bool test4b() {
    const double a_in=0.05, e_in=0.9, a_out=1.0;
    const double T_in = 2*M_PI*std::sqrt(a_in*a_in*a_in);
    // r_peri = a_in*(1-e_in) = 0.05*0.1 = 0.005
    // r_ks   = alpha_ks * a_semi ~ 0.6 * 0.05 = 0.03 > r_peri ✓
    ARChainNKSIntegrator::KSParameters ks;
    ks.alpha_ks          = 0.6;
    ks.hysteresis_factor = 2.0;
    ks.ks_substeps       = 64;
    ARChainNBSIntegrator::BSParameters bs;
    bs.bs_eps   = 1e-10;
    bs.max_steps = 500000;
    ARChainNKSIntegrator integ(1e-3, bs, ks);

    auto sys = make_triple_compact(0.5, 0.5, 0.2, a_in, e_in, a_out);
    const std::vector<int> idx = {0,1,2};
    const double E0 = compute_energy(sys);
    ARChainNState st = integ.initialize(sys, idx);

    const int    n  = 20;
    const double dt = T_in / n;
    bool activated  = false;

    auto t_start = std::chrono::steady_clock::now();

    for (int s = 0; s < n; ++s) {
        integ.integrate_to_ks(st, st.t_phys + dt);
        if (integ.ks_was_active()) activated = true;

        // Check de timeout en cada paso
        auto elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - t_start).count();
        if (elapsed > 60.0) {
            std::cout << "    [TIMEOUT] paso " << s+1 << "/" << n
                      << "  t_phys=" << st.t_phys
                      << "  elapsed=" << elapsed << "s\n";
            return FAIL("triple e=0.9 TIMEOUT",
                        "supero 60s — paso adaptativo KS insuficiente en triple "
                        "(ver Maestro v11 sec.7.3: fix en integrate_until())");
        }
    }

    auto elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_start).count();

    integ.write_back(st, sys, idx);
    const double dE = std::abs(compute_energy(sys)-E0)/std::abs(E0);

    std::cout << "    a_in=" << a_in << "  e_in=" << e_in
              << "  r_peri=" << (a_in*(1-e_in))
              << "  KS_act=" << (activated?"SI":"NO")
              << "  link_crit=" << integ.ks_critical_link()
              << "  |dE/E0|=" << dE
              << "  tiempo=" << std::fixed << std::setprecision(1)
              << elapsed << "s\n";

    if (!activated)
        return FAIL("triple e=0.9", "KS no se activo en el par interno");
    if (dE > 1e-1)
        return FAIL("triple e=0.9", "|dE|>1e-1 — divergencia");
    return PASS("triple (e=0.9): KS activo, sin divergencia, < 60s");
}

// T5 — Histeresis: transiciones ON/OFF <= 4 en 2 periodos (e=0.5)
static bool test5() {
    const double a=1.0,e=0.5,T=2*M_PI*std::sqrt(a*a*a);
    ARChainNKSIntegrator::KSParameters ks; ks.alpha_ks=0.6; ks.hysteresis_factor=2.0; ks.ks_substeps=64;
    ARChainNBSIntegrator::BSParameters bs; bs.bs_eps=1e-10; bs.max_steps=2000000;
    ARChainNKSIntegrator integ(1e-3,bs,ks);
    auto sys=make_binary_ecc(0.5,0.5,a,e);
    const std::vector<int> idx={0,1};
    ARChainNState st=integ.initialize(sys,idx);
    const int n=40; const double dt=T*2/n; int trans=0; bool prev=false;
    for(int s=0;s<n;++s){
        integ.integrate_to_ks(st,st.t_phys+dt);
        bool curr=integ.ks_was_active();
        if(curr!=prev)++trans; prev=curr;
    }
    std::cout<<"    e=0.5  alpha=0.6  transiciones="<<trans<<"  (2 periodos, 40 pasos)\n";
    if(trans>4) return FAIL("histeresis","transiciones>4 — flip-flopping detectado");
    return PASS("histeresis: transiciones <= 4 en 2 periodos, sin flip-flopping");
}

int main() {
    std::cout<<std::scientific<<std::setprecision(4);
    std::cout<<"======================================================\n";
    std::cout<<"  test_archain_n_ks — KS en AR-chain (Fase 6B)\n";
    std::cout<<"======================================================\n\n";
    int passed=0,total=0;
    auto run=[&](bool(*fn)(),const char* name){
        ++total;
        std::cout<<"Test "<<total<<" — "<<name<<"\n";
        try{if(fn())++passed;}catch(const std::exception& ex){FAIL(name,std::string("excepcion: ")+ex.what());}
        std::cout<<"\n";
    };
    run(test1,  "Circular (e=0): KS inactivo");
    run(test2,  "Excentrica (e=0.9): KS activo en perihelio");
    run(test3,  "Alta excentricidad (e=0.999): conservacion con KS");
    run(test4,  "Triple compacto (e=0.5): KS en par interno [base]");
    run(test4b, "Triple compacto (e=0.9): rendimiento < 60s [deuda tecnica]");
    run(test5,  "Histeresis: sin flip-flopping ON/OFF");
    std::cout<<"======================================================\n";
    std::cout<<"  Resultado: "<<passed<<"/"<<total<<" tests\n";
    std::cout<<"======================================================\n";
    return (passed==total)?0:1;
}