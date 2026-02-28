// core/state/State.h
// ============================================================================
// INTERFAZ CONCEPTUAL: State
// ============================================================================
// Un State representa el estado físico completo de un sistema en un instante t.
// Es datos puros — sin lógica de evolución temporal embebida.
//
// Principio de diseño (Guía Maestra):
//   - State contiene solo posiciones, velocidades y masas (o equivalentes).
//   - La evolución dState/dt es responsabilidad exclusiva del Integrator.
//   - Los observables (energía, momento angular) los calculan los Analyzers.
//
// Implementaciones concretas en este proyecto:
//   - Body            → estado de una partícula puntual (pos, vel, masa)
//   - NBodySystem     → colección de Body
//   - BinaryState     → estado de un par (posición/velocidad relativa + CM)
//   - KSState         → estado KS de una binaria (u, w, energy, tau)
//   - Chain3State     → estado Chain KS de un triple (Q1,P1,Q2,P2, CM, tau)
//
// No se define una clase base virtual State porque las implementaciones
// tienen tipos fundamentalmente distintos y la genericidad se logra mediante
// templates o sobrecarga de Integrator.
// ============================================================================
#pragma once
