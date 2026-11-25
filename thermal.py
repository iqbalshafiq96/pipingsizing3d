"""
thermal.py - Temperature drop along piping
==========================================
Heat loss from a (optionally insulated) pipe to its surroundings, plus the
temperature change caused by the pressure change itself:

  * Heat loss:   series thermal resistances per metre of pipe
                 R' = R_inside + R_wall + R_insulation + R_outside   [K.m/W]
                 T_out = T_amb + (T_in - T_amb) * exp(-dL / (R' * m * cp))
  * Gas:         Joule-Thomson effect      dT = mu_JT * dP
  * Liquid:      frictional heating        dT = dP_loss / (rho * cp)

All functions work in SI units (K, Pa, m, kg/s, W).
This module has no Streamlit / HTML dependency - hydraulics.py calls it.
"""
import math

SIGMA = 5.670374419e-8  # Stefan-Boltzmann, W/m2K4
G = 9.80665

# --------------------------------------------------------------------------
# Catalogues (sent to index.html so the drop-downs are data driven)
# --------------------------------------------------------------------------
INSULATION = {
    "None":                   {"k": None},
    "Mineral wool":           {"k": 0.040},
    "Glass wool":             {"k": 0.038},
    "Calcium silicate":       {"k": 0.055},
    "Cellular glass":         {"k": 0.045},
    "Polyurethane foam":      {"k": 0.025},
    "Polyisocyanurate (PIR)": {"k": 0.027},
    "Elastomeric foam":       {"k": 0.036},
    "Aerogel blanket":        {"k": 0.015},
    "Custom":                 {"k": None},
}

JACKETS = {  # outer surface emissivity
    "Aluminium jacket": 0.10,
    "Stainless jacket": 0.20,
    "Painted / PVC jacket": 0.90,
    "Bare (oxidised) steel": 0.80,
}

ENVIRONMENTS = {
    "air_still": "Air - still / indoor",
    "air_wind": "Air - outdoor with wind",
    "buried": "Buried in soil",
}


# --------------------------------------------------------------------------
# Inside film coefficient
# --------------------------------------------------------------------------
def inside_htc(Re, Pr, k_fluid, D, f_darcy):
    """Inside film coefficient h_i [W/m2K].
    Laminar Nu = 3.66, turbulent Gnielinski, linear blend 2300 < Re < 3000."""
    def gnielinski(re):
        f = f_darcy if f_darcy and f_darcy > 0 else (0.790 * math.log(re) - 1.64) ** -2
        return (f / 8) * (re - 1000) * Pr / (1 + 12.7 * math.sqrt(f / 8) * (Pr ** (2 / 3) - 1))

    if Re <= 2300:
        Nu = 3.66
    elif Re < 3000:
        w = (Re - 2300) / 700
        Nu = (1 - w) * 3.66 + w * max(gnielinski(3000), 3.66)
    else:
        Nu = max(gnielinski(Re), 3.66)
    return Nu * k_fluid / D


# --------------------------------------------------------------------------
# Outside coefficient (air) - natural/forced convection + radiation
# --------------------------------------------------------------------------
def _air_props(T_film):
    """Approximate dry-air properties at 1 atm as function of film temperature [K]."""
    k = 0.0241 * (T_film / 273.15) ** 0.81
    nu = 1.338e-5 * (T_film / 273.15) ** 1.75
    return k, nu, 0.71


def outside_htc_air(Ts, Ta, Do, wind, emissivity):
    """Combined convection + radiation coefficient for a cylinder in air."""
    Tf = 0.5 * (Ts + Ta)
    k, nu, Pr = _air_props(Tf)
    if wind and wind > 0.05:
        # Churchill-Bernstein (forced convection, cross flow)
        Re = wind * Do / nu
        Nu = 0.3 + (0.62 * Re ** 0.5 * Pr ** (1 / 3)) / (1 + (0.4 / Pr) ** (2 / 3)) ** 0.25 \
            * (1 + (Re / 282000) ** 0.625) ** 0.8
    else:
        # Churchill-Chu (natural convection, horizontal cylinder)
        dT = max(abs(Ts - Ta), 0.01)
        Ra = G * (1 / Tf) * dT * Do ** 3 / (nu * nu) * Pr
        Nu = (0.60 + 0.387 * Ra ** (1 / 6) / (1 + (0.559 / Pr) ** (9 / 16)) ** (8 / 27)) ** 2
    h_conv = Nu * k / Do
    h_rad = emissivity * SIGMA * (Ts ** 2 + Ta ** 2) * (Ts + Ta)
    return h_conv + h_rad, h_conv, h_rad


# --------------------------------------------------------------------------
# Overall resistance per metre
# --------------------------------------------------------------------------
def resistance_per_length(T_fluid, h_i, Di, Do, k_pipe, cfg):
    """Return dict with R' [K.m/W] and breakdown.
    cfg (SI): ins_k, ins_thk, env, T_amb, wind, emissivity, depth, k_soil"""
    R_in = 1.0 / (h_i * math.pi * Di) if h_i > 0 else 0.0
    R_wall = math.log(Do / Di) / (2 * math.pi * k_pipe) if Do > Di else 0.0
    ins_k, thk = cfg.get("ins_k"), cfg.get("ins_thk", 0.0)
    if ins_k and thk > 0:
        Dins = Do + 2 * thk
        R_ins = math.log(Dins / Do) / (2 * math.pi * ins_k)
    else:
        Dins, R_ins = Do, 0.0
    Ta = cfg["T_amb"]
    R_cond = R_in + R_wall + R_ins

    if cfg["env"] == "buried":
        H = max(cfg.get("depth", 1.0), 0.51 * Dins)
        R_out = math.acosh(2 * H / Dins) / (2 * math.pi * cfg.get("k_soil", 1.2))
        Ts, h_o = None, None
    else:
        wind = cfg.get("wind", 0.0) if cfg["env"] == "air_wind" else 0.0
        Ts = Ta + 0.3 * (T_fluid - Ta)  # first guess
        h_o = 10.0
        for _ in range(30):
            h_o, _, _ = outside_htc_air(Ts, Ta, Dins, wind, cfg.get("emissivity", 0.9))
            R_out = 1.0 / (h_o * math.pi * Dins)
            q = (T_fluid - Ta) / (R_cond + R_out)
            Ts_new = T_fluid - q * R_cond
            if abs(Ts_new - Ts) < 0.01:
                Ts = Ts_new
                break
            Ts = 0.5 * Ts + 0.5 * Ts_new
        R_out = 1.0 / (h_o * math.pi * Dins)

    R_tot = R_cond + R_out
    return {"R": R_tot, "R_in": R_in, "R_wall": R_wall, "R_ins": R_ins, "R_out": R_out,
            "Ts": Ts, "h_o": h_o, "Dins": Dins, "U_od": 1.0 / (R_tot * math.pi * Do)}


def heat_loss_step(T_in, m, cp, R_per_m, dL, T_amb):
    """Exponential (log-mean) heat-loss step. Returns (T_out, Q_lost[W])."""
    if m <= 0 or cp <= 0 or R_per_m <= 0:
        return T_in, 0.0
    T_out = T_amb + (T_in - T_amb) * math.exp(-dL / (R_per_m * m * cp))
    return T_out, m * cp * (T_in - T_out)


def pressure_effect_dT(phase, dP_total, dP_loss, rho, cp, mu_jt):
    """Temperature change caused by the pressure change.
    dP_total = P_out - P_in (negative for a drop), dP_loss = irreversible loss (positive).
    Gas   : isenthalpic Joule-Thomson   dT = mu_JT * dP_total
    Liquid: dissipated friction -> heat  dT = dP_loss / (rho*cp)"""
    if phase == "gas":
        return mu_jt * dP_total
    return dP_loss / (rho * cp) if rho > 0 and cp > 0 else 0.0


def thermal_config(th, T_amb_K):
    """Convert the thermal block from index.html into an SI config dict."""
    ins = th.get("insulation", "None")
    ins_k = th.get("k_custom") if ins == "Custom" else INSULATION.get(ins, {}).get("k")
    thk = float(th.get("thk_mm", 0) or 0) / 1000.0
    if not ins_k or ins == "None":
        ins_k, thk = None, 0.0
    if ins_k is None or thk == 0:
        emis = JACKETS.get("Bare (oxidised) steel")
    else:
        emis = JACKETS.get(th.get("jacket", "Aluminium jacket"), 0.1)
    return {
        "enabled": bool(th.get("enabled", True)),
        "ins_k": float(ins_k) if ins_k else None,
        "ins_thk": thk,
        "env": th.get("env", "air_still"),
        "T_amb": T_amb_K,
        "wind": float(th.get("wind", 0) or 0),
        "emissivity": emis,
        "depth": float(th.get("depth", 1.0) or 1.0),
        "k_soil": float(th.get("k_soil", 1.2) or 1.2),
    }


def catalog():
    return {"insulation": {k: v["k"] for k, v in INSULATION.items()},
            "jackets": JACKETS, "environments": ENVIRONMENTS}
