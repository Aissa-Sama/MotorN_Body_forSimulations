// python/nbody_bindings.cpp
// ============================================================================
// PYBIND11 BINDINGS — nbody_core — Fases 1–7D
//
// Cambios respecto a la version 6D:
//   - step_to(): eliminado dt_hint (firma actual: system, t_final, on_step)
//   - pn_active_count(): eliminado (no existe en HierarchicalIntegrator)
//   - BodyType: expuesto como enum (Fase 7C)
//   - Body.radius, Body.type: expuestos (Fase 7C)
//   - Trajectory + integrate(): arrays NumPy (n_snapshots, N, 3) (Fase 7D)
// ============================================================================
#define _USE_MATH_DEFINES
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <cmath>
#include <vector>

#include "vec3.h"
#include "body.h"
#include "nbody_system.h"
#include "initial_conditions.h"
#include "leapfrog_integrator.h"
#include "hierarchical_integrator.h"
#include "hierarchy_builder.h"
#include "physics_observables.h"

namespace py = pybind11;
using namespace pybind11::literals;

// ============================================================================
// TRAJECTORY — Fase 7D
// ============================================================================
struct Trajectory {
    py::array_t<double> positions;   // (n_snapshots, N, 3)
    py::array_t<double> velocities;  // (n_snapshots, N, 3)
    py::array_t<double> times;       // (n_snapshots,)
    py::array_t<double> energies;    // (n_snapshots,)
};

static Trajectory integrate_trajectory(
    HierarchicalIntegrator& integ,
    NBodySystem& sys,
    double t_final,
    int n_snapshots,
    double dt_step)
{
    const int N = static_cast<int>(sys.bodies.size());
    const double dt_snap = t_final / static_cast<double>(n_snapshots);

    std::vector<double> pos_data, vel_data, times_data, energies_data;
    pos_data.reserve(n_snapshots * N * 3);
    vel_data.reserve(n_snapshots * N * 3);
    times_data.reserve(n_snapshots);
    energies_data.reserve(n_snapshots);

    std::vector<bool> used(N, false);
    double t_cur = 0.0;

    for (int s = 0; s < n_snapshots; ++s) {
        const double t_next = (s + 1) * dt_snap;
        while (t_cur < t_next - 1e-14) {
            const double dt = std::min(dt_step, t_next - t_cur);
            integ.step(sys, dt, used);
            t_cur += dt;
        }
        times_data.push_back(t_cur);
        energies_data.push_back(sys.total_energy());
        for (const auto& b : sys.bodies) {
            pos_data.push_back(b.position.x);
            pos_data.push_back(b.position.y);
            pos_data.push_back(b.position.z);
            vel_data.push_back(b.velocity.x);
            vel_data.push_back(b.velocity.y);
            vel_data.push_back(b.velocity.z);
        }
    }

    Trajectory traj;
    traj.positions  = py::array_t<double>(
        {n_snapshots, N, 3}, pos_data.data());
    traj.velocities = py::array_t<double>(
        {n_snapshots, N, 3}, vel_data.data());
    traj.times      = py::array_t<double>(
        {n_snapshots},       times_data.data());
    traj.energies   = py::array_t<double>(
        {n_snapshots},       energies_data.data());
    return traj;
}

PYBIND11_MODULE(nbody_core, m) {
    m.doc() = "nbody_core — Motor N-body gravitacional de alta precision (C++17)";

    // ========================================================================
    // Vec3
    // ========================================================================
    py::class_<Vec3>(m, "Vec3")
        .def(py::init<>())
        .def(py::init<double, double, double>(), py::arg("x"), py::arg("y"), py::arg("z"))
        .def_readwrite("x", &Vec3::x)
        .def_readwrite("y", &Vec3::y)
        .def_readwrite("z", &Vec3::z)
        .def("norm",  &Vec3::norm,  "Norma euclidiana |v|")
        .def("norm2", &Vec3::norm2, "Norma al cuadrado |v|^2")
        .def("dot",   &Vec3::dot,   "Producto punto", py::arg("other"))
        .def(py::self + py::self)
        .def(py::self - py::self)
        .def(py::self * double())
        .def(double() * py::self)
        .def(py::self / double())
        .def(-py::self)
        .def("__repr__", [](const Vec3& v) {
            return "Vec3(" + std::to_string(v.x) + ", "
                           + std::to_string(v.y) + ", "
                           + std::to_string(v.z) + ")";
        });

    m.def("dot",   &dot,   "Producto punto de dos Vec3",    py::arg("a"), py::arg("b"));
    m.def("cross", &cross, "Producto vectorial de dos Vec3", py::arg("a"), py::arg("b"));

    // ========================================================================
    // BodyType (Fase 7C)
    // ========================================================================
    py::enum_<BodyType>(m, "BodyType")
        .value("POINT_MASS",     BodyType::POINT_MASS)
        .value("STAR",           BodyType::STAR)
        .value("COMPACT_OBJECT", BodyType::COMPACT_OBJECT)
        .value("PLANET",         BodyType::PLANET)
        .value("SATELLITE",      BodyType::SATELLITE)
        .export_values();

    // ========================================================================
    // Body
    // ========================================================================
    py::class_<Body>(m, "Body")
        .def(py::init<>())
        .def_readwrite("position", &Body::position)
        .def_readwrite("velocity", &Body::velocity)
        .def_readwrite("mass",     &Body::mass)
        .def_readwrite("radius",   &Body::radius)
        .def_readwrite("type",     &Body::type)
        .def("is_extended",   &Body::is_extended)
        .def("is_compact",    &Body::is_compact)
        .def("is_star",       &Body::is_star)
        .def("is_point_mass", &Body::is_point_mass)
        .def("__repr__", [](const Body& b) {
            return "Body(m=" + std::to_string(b.mass)
                + ", r=(" + std::to_string(b.position.x) + ", "
                + std::to_string(b.position.y) + ", "
                + std::to_string(b.position.z) + ")"
                + ", radius=" + std::to_string(b.radius) + ")";
        });

    // ========================================================================
    // NBodySystem
    // ========================================================================
    py::class_<NBodySystem>(m, "NBodySystem")
        .def(py::init<>())
        .def_readwrite("bodies", &NBodySystem::bodies)
        .def_readwrite("G",      &NBodySystem::G)
        .def("compute_accelerations",    &NBodySystem::compute_accelerations)
        .def("invalidate_accelerations", &NBodySystem::invalidate_accelerations)
        .def("kinetic_energy",           &NBodySystem::kinetic_energy)
        .def("potential_energy",         &NBodySystem::potential_energy)
        .def("total_energy",             &NBodySystem::total_energy)
        .def("total_momentum",           &NBodySystem::total_momentum)
        .def("total_angular_momentum",   &NBodySystem::total_angular_momentum)
        .def("__repr__", [](const NBodySystem& sys) {
            return "NBodySystem(N=" + std::to_string(sys.bodies.size())
                + ", G=" + std::to_string(sys.G) + ")";
        });

    // ========================================================================
    // InitialConditions
    // ========================================================================
    py::class_<InitialConditions>(m, "InitialConditions")
        .def_static("solar_system",      &InitialConditions::solar_system)
        .def_static("binary_with_field", &InitialConditions::binary_with_field)
        .def_static("plummer_cluster",   &InitialConditions::plummer_cluster,
                    py::arg("n_bodies"))
        .def_static("figure_eight",      &InitialConditions::figure_eight)
        .def_static("random_system",     &InitialConditions::random_system,
                    py::arg("n_bodies"), py::arg("radius"))
        .def_static("kepler_binary",     &InitialConditions::kepler_binary,
                    py::arg("a") = 1.0, py::arg("e") = 0.0);

    // ========================================================================
    // HierarchyBuilder::Params
    // ========================================================================
    py::class_<HierarchyBuilder::Params::PNParams>(m, "PNParams")
        .def(py::init<>())
        .def_readwrite("enabled",            &HierarchyBuilder::Params::PNParams::enabled)
        .def_readwrite("c_speed",            &HierarchyBuilder::Params::PNParams::c_speed)
        .def_readwrite("pn_order",           &HierarchyBuilder::Params::PNParams::pn_order)
        .def_readwrite("activation_epsilon", &HierarchyBuilder::Params::PNParams::activation_epsilon);

    py::class_<HierarchyBuilder::Params>(m, "HierarchyParams")
        .def(py::init<>())
        .def_readwrite("r_ks_threshold",      &HierarchyBuilder::Params::r_ks_threshold)
        .def_readwrite("tidal_threshold",     &HierarchyBuilder::Params::tidal_threshold)
        .def_readwrite("strong_coupling_eta", &HierarchyBuilder::Params::strong_coupling_eta)
        .def_readwrite("ar_chain_threshold",  &HierarchyBuilder::Params::ar_chain_threshold)
        .def_readwrite("ar_chain_eta",        &HierarchyBuilder::Params::ar_chain_eta)
        .def_readwrite("pn",                  &HierarchyBuilder::Params::pn);

    // ========================================================================
    // HierarchicalIntegrator
    // ========================================================================
    py::class_<HierarchicalIntegrator>(m, "HierarchicalIntegrator")
        .def(py::init([](double r_ks, double ks_dt,
                         const HierarchyBuilder::Params& params) {
            return std::make_unique<HierarchicalIntegrator>(
                std::make_unique<LeapfrogIntegrator>(),
                r_ks, ks_dt, params, nullptr
            );
        }),
            py::arg("r_ks_threshold") = 1.0,
            py::arg("ks_internal_dt") = 1e-4,
            py::arg("params") = HierarchyBuilder::Params{}
        )
        .def("step", [](HierarchicalIntegrator& self,
                        NBodySystem& system, double dt) {
            std::vector<bool> used(system.bodies.size(), false);
            self.step(system, dt, used);
        }, py::arg("system"), py::arg("dt"))
        .def("step_to",
            [](HierarchicalIntegrator& self,
               NBodySystem& sys, double t_final,
               std::function<void(double)> cb) {
                self.step_to(sys, t_final, cb);
            },
            py::arg("system"),
            py::arg("t_final"),
            py::arg("on_step") = std::function<void(double)>{})
        // integrate: Fase 7D — devuelve trayectoria completa como numpy arrays
        .def("integrate", &integrate_trajectory,
             py::arg("system"),
             py::arg("t_final"),
             py::arg("n_snapshots") = 100,
             py::arg("dt_step")     = 0.01,
             R"pbdoc(
Integra hasta t_final y devuelve una Trajectory con arrays NumPy.

Parametros
----------
system      : NBodySystem  (modificado in-place)
t_final     : float        tiempo fisico final
n_snapshots : int          numero de snapshots uniformemente espaciados
dt_step     : float        paso maximo de integracion entre snapshots

Retorna
-------
Trajectory con:
  .positions   ndarray (n_snapshots, N, 3)
  .velocities  ndarray (n_snapshots, N, 3)
  .times       ndarray (n_snapshots,)
  .energies    ndarray (n_snapshots,)
)pbdoc");

    // ========================================================================
    // Trajectory (Fase 7D)
    // ========================================================================
    py::class_<Trajectory>(m, "Trajectory")
        .def_readwrite("positions",  &Trajectory::positions)
        .def_readwrite("velocities", &Trajectory::velocities)
        .def_readwrite("times",      &Trajectory::times)
        .def_readwrite("energies",   &Trajectory::energies)
        .def("__repr__", [](const Trajectory& t) {
            return "<Trajectory snapshots=" +
                   std::to_string(t.times.size()) + ">";
        })
        .def("n_snapshots", [](const Trajectory& t) {
            return t.times.size();
        })
        .def("n_bodies", [](const Trajectory& t) {
            auto info = t.positions.request();
            return info.ndim >= 2 ? info.shape[1] : 0;
        });

    // ========================================================================
    // Physics Observables (namespace phys)
    // ========================================================================
    py::class_<phys::Observables>(m, "Observables")
        .def(py::init<>())
        .def_readonly("t",           &phys::Observables::t)
        .def_readonly("E",           &phys::Observables::E)
        .def_readonly("E_kinetic",   &phys::Observables::E_kinetic)
        .def_readonly("E_potential", &phys::Observables::E_potential)
        .def_readonly("dE_rel",      &phys::Observables::dE_rel)
        .def_readonly("P",           &phys::Observables::P)
        .def_readonly("dP_rel",      &phys::Observables::dP_rel)
        .def_readonly("L",           &phys::Observables::L)
        .def_readonly("dL_rel",      &phys::Observables::dL_rel)
        .def_readonly("cm_pos",      &phys::Observables::cm_pos)
        .def_readonly("cm_vel",      &phys::Observables::cm_vel)
        .def_readonly("sep_min",     &phys::Observables::sep_min)
        .def_readonly("sep_max",     &phys::Observables::sep_max)
        .def_readonly("virial",      &phys::Observables::virial)
        .def_readonly("regime",      &phys::Observables::regime);

    m.def("kinetic_energy",         &phys::kinetic_energy,         py::arg("system"));
    m.def("potential_energy",       &phys::potential_energy,
          py::arg("system"), py::arg("softening") = 0.0);
    m.def("total_energy",           &phys::total_energy,
          py::arg("system"), py::arg("softening") = 0.0);
    m.def("total_momentum",         &phys::total_momentum,         py::arg("system"));
    m.def("total_angular_momentum", &phys::total_angular_momentum, py::arg("system"));
    m.def("total_mass",             &phys::total_mass,             py::arg("system"));
    m.def("min_separation",         &phys::min_separation,         py::arg("system"));
    m.def("max_separation",         &phys::max_separation,         py::arg("system"));
    m.def("virial_ratio",           &phys::virial_ratio,           py::arg("system"));
    m.def("compute_observables",    &phys::compute_observables,
          py::arg("system"), py::arg("t"), py::arg("E0"),
          py::arg("P0"), py::arg("L0"),
          py::arg("regime") = "", py::arg("softening") = 0.0);
}
