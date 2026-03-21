#pragma once
// ============================================================================
// analysis/merger_engine.h
//
// Ejecutor de fusiones y eventos de disruption tidal.
// Opera sobre NBodySystem directamente, modificando bodies[].
//
// CONSERVACIÓN:
//   - Colisión física:    masa, momento lineal, posición CM — exactos.
//   - Full TDE:           BH gana fracción f_acc de la masa estelar.
//                         La fracción (1-f_acc) sale del sistema como "debris".
//                         Masa, momentum del BH actualizados exactamente.
//   - Partial TDE:        estrella pierde fracción de masa proporcional a
//                         la penetración en r_t (modelo de Stone et al. 2013).
//   - Schwarzschild cap:  cuerpo ligero acretado completamente por el BH.
//
// MODELO DE ACRECIÓN TDE (Rees 1988; Evans & Kochanek 1989):
//   f_acc = 0.5 exacto para TDE clásico en parábola.
//   El debris que escapa lleva el momento del cuerpo ligero menos f_acc.
//   En el motor: el debris se elimina del sistema (o se registra como
//   un cuerpo fantasma de masa residual, según simplified_mergers).
//
// PARTIAL TDE — modelo de Stone, Sari & Loeb (2013):
//   dm/m = (1 - (r_p/r_t)^(3/2))    para r_t/2 < r_p < r_t
//   El radio de la estrella se reduce: R_new = R_old * (m_new/m_old)^(1/3)
//   (conserva densidad media)
//
// FLAG simplified_mergers:
//   true  → TDE = fusión instantánea (BH absorbe todo). Para pruebas rápidas.
//   false → física completa (default).
//
// REFERENCIAS:
//   Rees (1988), Nature 333, 523
//   Evans & Kochanek (1989), ApJ 346, L13
//   Stone, Sari & Loeb (2013), MNRAS 435, 1809
// ============================================================================
#include "body.h"
#include "tidal_events.h"
#include "nbody_system.h"
#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include <set>

namespace tidal {

// ── Resultado de una fusión ───────────────────────────────────────────────────

struct MergerResult {
    bool   occurred      = false;
    int    survivor_idx  = -1;   // índice del cuerpo remanente en bodies[]
    int    removed_idx   = -1;   // índice eliminado
    double mass_accreted = 0.0;  // masa ganada por el survivor
    double mass_lost     = 0.0;  // masa que salió del sistema (debris TDE)
    double delta_L_orb   = 0.0;  // momento angular orbital depositado (J)
    EventType type       = EventType::NONE;
};

// ── Callback de logging ───────────────────────────────────────────────────────
// Llamado cada vez que ocurre un evento. Firma:
//   void on_event(const TidalEvent&, const MergerResult&)
using MergerCallback = std::function<void(const TidalEvent&, const MergerResult&)>;

// ── MergerEngine ──────────────────────────────────────────────────────────────

class MergerEngine {
public:

    struct Params {
        double f_accretion     = 0.5;    // fracción de masa estelar acretada en Full TDE
        double c_light         = 0.0;    // velocidad de la luz en N-body units (0=no relativista)
        bool   simplified      = false;  // true → TDE = fusión instantánea
        bool   log_events      = true;   // imprimir eventos a stdout
        bool   keep_debris     = false;  // true → debris TDE permanece como cuerpo de baja masa
        double debris_min_mass = 1e-6;   // masa mínima para mantener debris como cuerpo
    };

    explicit MergerEngine(const Params& p = Params{}) : params_(p) {}

    void set_callback(MergerCallback cb) { callback_ = cb; }

    // ── Interfaz principal ────────────────────────────────────────────────────
    //
    // Detecta y ejecuta todos los eventos presentes en el sistema.
    // Modifica system.bodies[] in-place.
    // Devuelve número de eventos ejecutados.
    //
    int process(NBodySystem& system, double t_now = 0.0) {
        int n_events = 0;
        bool any = true;

        // Pares que ya recibieron un PARTIAL_TDE en esta llamada.
        // No se redetectan dentro del mismo paso — el par necesita
        // separarse orbitalmente antes de poder ser afectado de nuevo.
        std::set<std::pair<int,int>> partial_tde_done;

        while (any) {
            auto events = tidal::detect_events(system.bodies, params_.c_light, t_now);
            any = false;

            for (const auto& ev : events) {
                if (ev.i >= (int)system.bodies.size() ||
                    ev.j >= (int)system.bodies.size()) continue;

                // No redetectar PARTIAL_TDE sobre el mismo par en este paso
                if (ev.type == tidal::EventType::PARTIAL_TDE) {
                    auto key = std::make_pair(std::min(ev.i,ev.j),
                                             std::max(ev.i,ev.j));
                    if (partial_tde_done.count(key)) continue;
                    partial_tde_done.insert(key);
                }

                MergerResult result = execute(system, ev);
                if (result.occurred) {
                    ++n_events;
                    // Solo reiniciar el bucle si hubo eliminación de cuerpo
                    // (colisión, full TDE, captura) — no para partial TDE
                    if (result.removed_idx >= 0) any = true;
                    if (callback_) callback_(ev, result);
                    if (params_.log_events) log_event(ev, result, t_now);
                    if (result.removed_idx >= 0) break; // redetectar con índices actualizados
                }
            }
        }

        return n_events;
    }

private:
    Params          params_;
    MergerCallback  callback_;

    // ── Ejecutar un evento individual ─────────────────────────────────────────

    MergerResult execute(NBodySystem& sys, const TidalEvent& ev) {
        switch (ev.type) {
            case EventType::PHYSICAL_COLLISION:
                return execute_collision(sys, ev.i, ev.j);
            case EventType::FULL_TDE:
                return execute_full_tde(sys, ev.i, ev.j);
            case EventType::PARTIAL_TDE:
                if (params_.simplified)
                    return execute_full_tde(sys, ev.i, ev.j);  // simplificado
                return execute_partial_tde(sys, ev.i, ev.j, ev.sep, ev.r_crit);
            case EventType::SCHWARZSCHILD_CAP:
                return execute_schwarzschild_cap(sys, ev.i, ev.j);
            default:
                return MergerResult{};
        }
    }

    // ── Colisión física: fusión conservativa ──────────────────────────────────
    //
    // Conserva exactamente: masa total, momento lineal total, posición CM.
    // El momento angular orbital se deposita como spin del remanente
    // (registrado en delta_L_orb pero no modelado internamente).
    //
    MergerResult execute_collision(NBodySystem& sys, int i, int j) {
        auto& a = sys.bodies[i];
        auto& b = sys.bodies[j];

        const double M = a.mass + b.mass;

        // Momento angular orbital del par en el frame del CM
        Vec3 r_rel = b.position - a.position;
        Vec3 v_rel = b.velocity - a.velocity;
        double mu  = a.mass * b.mass / M;
        Vec3   L   = cross(r_rel, v_rel) * mu;
        double L_norm = L.norm();

        // Posición y velocidad del remanente — CM exacto
        Vec3 R_new = (a.position * a.mass + b.position * b.mass) * (1.0 / M);
        Vec3 V_new = (a.velocity * a.mass + b.velocity * b.mass) * (1.0 / M);

        // Radio del remanente: conserva volumen (densidad media constante)
        double R_phys = 0.0;
        if (a.is_extended() && b.is_extended())
            R_phys = std::cbrt(a.radius*a.radius*a.radius + b.radius*b.radius*b.radius);
        else if (a.is_extended()) R_phys = a.radius;
        else if (b.is_extended()) R_phys = b.radius;

        // Tipo del remanente: el más "compacto" gana
        BodyType T_new = (static_cast<int>(a.type) > static_cast<int>(b.type)) ? a.type : b.type;

        // Remanente en posición i, eliminar j
        a.mass     = M;
        a.position = R_new;
        a.velocity = V_new;
        a.radius   = R_phys;
        a.type     = T_new;

        remove_body(sys, j);

        MergerResult res;
        res.occurred      = true;
        res.survivor_idx  = i;
        res.removed_idx   = j;
        res.mass_accreted = b.mass;  // antes de eliminar
        res.mass_lost     = 0.0;
        res.delta_L_orb   = L_norm;
        res.type          = EventType::PHYSICAL_COLLISION;
        return res;
    }

    // ── Full TDE: BH disrumpe estrella ────────────────────────────────────────
    //
    // BH acrece f_acc de la masa estelar.
    // (1 - f_acc) sale como debris (eliminado o mantenido como cuerpo ligero).
    //
    // El momentum del BH se actualiza conservando el momentum total:
    //   p_BH_new = p_BH + p_star (todo el momentum de la estrella va al BH/debris)
    //   Dado que debris lleva (1-f_acc)*p_star y BH lleva p_BH + f_acc*p_star:
    //   V_BH_new = (M_BH*V_BH + f_acc*m_star*V_star) / (M_BH + f_acc*m_star)
    //
    MergerResult execute_full_tde(NBodySystem& sys, int i, int j) {
        // Identificar quién es el BH y quién la estrella
        int idx_bh   = (sys.bodies[i].mass >= sys.bodies[j].mass) ? i : j;
        int idx_star = (sys.bodies[i].mass >= sys.bodies[j].mass) ? j : i;

        auto& bh   = sys.bodies[idx_bh];
        auto& star = sys.bodies[idx_star];

        const double f   = params_.f_accretion;  // 0.5 por defecto
        const double dm  = f * star.mass;         // masa acretada
        const double m_debris = (1.0 - f) * star.mass;

        // Actualizar BH: masa y velocidad
        Vec3 V_bh_new = (bh.velocity * bh.mass + star.velocity * dm)
                       * (1.0 / (bh.mass + dm));
        bh.mass     += dm;
        bh.velocity  = V_bh_new;
        // Posición del BH no cambia (ya es el objeto masivo)

        MergerResult res;
        res.occurred      = true;
        res.survivor_idx  = idx_bh;
        res.removed_idx   = idx_star;
        res.mass_accreted = dm;
        res.mass_lost     = m_debris;
        res.type          = EventType::FULL_TDE;

        if (params_.keep_debris && m_debris > params_.debris_min_mass) {
            // Debris: 50% de la masa estelar, viaja en dirección del momento estelar residual
            auto& s = sys.bodies[idx_star];
            s.mass     = m_debris;
            s.type     = BodyType::POINT_MASS;
            s.radius   = 0.0;
            // velocidad del debris: conservación de momentum del sistema estrella
            // V_debris = (m_star*V_star - dm*V_star) / m_debris = V_star
            // (el debris viaja con la velocidad original de la estrella)
        } else {
            remove_body(sys, idx_star);
        }

        return res;
    }

    // ── Partial TDE: estrella pierde masa pero sobrevive ─────────────────────
    //
    // Régimen: r_t ≤ sep < 2·r_t  →  β = r_t/sep ∈ (0.5, 1.0)
    //
    // Guillochon & Ramirez-Ruiz (2013), ApJ 767, 25 — ec. de fracción de masa:
    //   f_loss ≈ max(0, β - β_min) / (1 - β_min)   con β_min ≈ 0.5
    //
    // Para β muy cercano a 1 (sep → r_t), f_loss → 1 y la estrella
    // pierde casi toda su masa (transición al Full TDE).
    // Para β = 0.5 (sep = 2·r_t), f_loss = 0 (sin pérdida de masa).
    //
    // NOTA: Stone, Sari & Loeb (2013) usan dm/m = 1-(1/β)^(3/2) para β>1,
    // que es el régimen Full TDE. Para β<1 (partial), la fórmula correcta
    // es la de Guillochon & Ramirez-Ruiz (2013).
    //
    MergerResult execute_partial_tde(NBodySystem& sys, int i, int j,
                                     double sep, double r_t) {
        int idx_bh   = (sys.bodies[i].mass >= sys.bodies[j].mass) ? i : j;
        int idx_star = (sys.bodies[i].mass >= sys.bodies[j].mass) ? j : i;

        auto& bh   = sys.bodies[idx_bh];
        auto& star = sys.bodies[idx_star];

        // β = r_t/sep ∈ (0.5, 1) para partial TDE
        const double beta    = r_t / sep;
        const double beta_min = 0.5;
        const double f_loss  = std::max(0.0, (beta - beta_min) / (1.0 - beta_min));
        const double dm_lost = f_loss * star.mass;
        const double dm_acc  = params_.f_accretion * dm_lost;

        // Actualizar estrella: pierde masa, radio se reduce (densidad constante)
        double m_star_new = star.mass - dm_lost;
        star.radius = star.radius * std::cbrt(m_star_new / star.mass);
        star.mass   = m_star_new;

        // Actualizar BH: acrece dm_acc
        Vec3 V_bh_new = (bh.velocity * bh.mass + star.velocity * dm_acc)
                       * (1.0 / (bh.mass + dm_acc));
        bh.mass     += dm_acc;
        bh.velocity  = V_bh_new;

        MergerResult res;
        res.occurred      = true;
        res.survivor_idx  = idx_star;  // estrella sobrevive (modificada)
        res.removed_idx   = -1;        // nadie eliminado
        res.mass_accreted = dm_acc;
        res.mass_lost     = dm_lost - dm_acc;
        res.type          = EventType::PARTIAL_TDE;
        return res;
    }

    // ── Captura relativista: cuerpo cae dentro del ISCO ───────────────────────
    //
    // El objeto ligero es acretado completamente por el BH.
    // Conserva masa y momento lineal exactamente.
    //
    MergerResult execute_schwarzschild_cap(NBodySystem& sys, int i, int j) {
        int idx_bh  = (sys.bodies[i].mass >= sys.bodies[j].mass) ? i : j;
        int idx_obj = (sys.bodies[i].mass >= sys.bodies[j].mass) ? j : i;

        auto& bh  = sys.bodies[idx_bh];
        auto& obj = sys.bodies[idx_obj];

        const double M = bh.mass + obj.mass;
        Vec3 V_new = (bh.velocity * bh.mass + obj.velocity * obj.mass) * (1.0 / M);

        bh.mass     = M;
        bh.velocity = V_new;

        remove_body(sys, idx_obj);

        MergerResult res;
        res.occurred      = true;
        res.survivor_idx  = idx_bh;
        res.removed_idx   = idx_obj;
        res.mass_accreted = obj.mass;
        res.mass_lost     = 0.0;
        res.type          = EventType::SCHWARZSCHILD_CAP;
        return res;
    }

    // ── Eliminar cuerpo del sistema ───────────────────────────────────────────
    //
    // Elimina bodies[idx] con swap-and-pop para evitar O(N) shifts.
    // Nota: esto cambia los índices de los cuerpos posteriores.
    // El bucle de process() redetecta desde cero después de cada evento.
    //
    void remove_body(NBodySystem& sys, int idx) {
        if (idx < 0 || idx >= (int)sys.bodies.size()) return;
        // swap con el último y pop
        sys.bodies[idx] = sys.bodies.back();
        sys.bodies.pop_back();
        sys.invalidate_accelerations();
    }

    // ── Log ───────────────────────────────────────────────────────────────────

    void log_event(const TidalEvent& ev, const MergerResult& res, double t) const {
        const char* type_str[] = { "NONE", "COLLISION", "FULL_TDE",
                                   "PARTIAL_TDE", "ROCHE_OVERFLOW", "SCHWARZSCHILD_CAP" };
        int ti = static_cast<int>(ev.type);
        std::cout << "[MERGER t=" << t << "] "
                  << type_str[ti]
                  << " i=" << ev.i << " j=" << ev.j
                  << " sep=" << ev.sep << " r_crit=" << ev.r_crit
                  << " dm_acc=" << res.mass_accreted
                  << " dm_lost=" << res.mass_lost << "\n";
    }
};

} // namespace tidal