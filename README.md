# nbody_core

A modular C++ N-body simulation core with hybrid integration and exact regularization for close encounters.

This project is focused on **physical correctness, long-term stability, and extensibility**, rather than real-time performance or visualization alone.

---

## Features

- General N-body gravitational dynamics
- Multiple integrators:
  - Euler
  - Leapfrog
  - RK4
  - Velocity Verlet
- Hybrid integration framework
- Exact Kustaanheimo–Stiefel (KS) regularization for bound binaries
- Regime-aware integration (global ↔ binary)
- Energy, momentum and angular momentum diagnostics
- Decoupled logging and testing infrastructure

---

## Project Status

The core physical architecture is **functional and validated**.

Current capabilities:
- Stable long-term integration of binaries using KS
- Hybrid integration without double-stepping
- Verified energy and phase conservation for regularized orbits

The project is under active development.

---

## Roadmap (Next Steps)

### 1. Binary isolation and conflict resolution
Implement robust detection and isolation of multiple bound binaries:
- Candidate list construction (`BinaryPair`)
- Greedy but consistent selection strategy
- Avoidance of overlapping binary assignments

### 2. Strict temporal consistency
Ensure:
- KS binaries advance exactly the same physical Δt as the global system
- No discontinuities when entering or exiting regularized regimes

### 3. Snapshot layer for visualization
Separate:
- Internal integration state (KS variables, internal coordinates)
- Observable physical state (positions, velocities)

This enables:
- Smooth visualization
- Clean data export
- Python bindings without coupling to integrators

### 4. Python bindings
Expose the simulation core for:
- Analysis
- Visualization
- Experimentation
- Educational use

---

## Philosophy

This project prioritizes:
- Physical realism
- Clear separation of concerns
- Testable numerical behavior
- Extensibility over premature optimization

It is designed as a **scientific core**, not a monolithic application.

---

## Build

Basic example (Windows / PowerShell):

```powershell
mkdir build
cd build
cmake ..
cmake --build .



The project includes numerical validation tests, including:
- KS energy conservation
- Orbital phase consistency

Example:
.\test_ks.exe


```text
nbody_core/
├── CMakeLists.txt
├── README.md
├── LICENSE
├── src/ (src/ no implemented)      
│   ├── analysis/
│   ├── build/
│   ├── out/
│   ├── core/
│   ├── integrators/
│   ├── regularization/
│   │   └── ks/
│   └── io/
├── tests/
│   └── test_ks.cpp
└── examples/   =====>> No yet