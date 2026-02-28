// regularization/ks/ks_perturbed_integrator.h
#pragma once
#include "ks_integrator.h"
#include "nbody_system.h"
#include "vec3.h"
#include <vector>

class KSPerturbedIntegrator : public KSIntegrator {
public:
    // Declaración explícita del constructor
    explicit KSPerturbedIntegrator(double internal_dt);

    void integrate_perturbed(
        BinaryState& binary,
        double dt_global,
        const NBodySystem& system,
        int idx_i, int idx_j
    );

private:
    struct Perturber {
        int index;
        double mass;
        Vec3 position;
        Vec3 velocity;
        double distance;
    };

    std::vector<Perturber> find_perturbers(
        const NBodySystem& system,
        int idx_i, int idx_j,
        double r_binary
    ) const;

    Vec3 compute_perturbation_force(
        const BinaryState& binary,
        const std::vector<Perturber>& perturbers,
        const NBodySystem& system
    ) const;

    void force_to_ks(const Vec3& F_phys, const double u[4], double F_ks[4]) const;

    void leapfrog_step_perturbed(
        KSState& ks,
        const double F_ks[4],
        double dtau
    ) const;
};