"""
hydraulics.py - Piping hydraulics engine (3D routing)
=====================================================
Takes the network (tree) built in index.html and returns, for every segment:

  Continuity      v  = m / (rho * A)
  Reynolds        Re = rho * v * D / mu         -> flow regime
  Friction factor Darcy f (laminar 64/Re, Colebrook-White turbulent, blended 2300-4000)
  Darcy-Weisbach  dP_f = f (L/D) rho v^2 / 2
  Minor losses    dP_k = K rho v^2 / 2          (Crane K = n * f_T)
  Energy balance  P_out = P_in - dP_f - dP_k - rho g dz - (rho/2)(v2^2 - v1^2)

Geometry is 3D: every segment runs East/West/North/South/Up/Down; horizontal pipes
may carry an elevation change dz (sloped line). Coordinates are reported in metres
(E = +x, N = +y, Z = elevation).
Temperature along the line is delegated to thermal.py.
Everything inside is SI; conversion happens only at the boundary.
"""
import copy
import math

import thermal

G = 9.80665
R_UNIV = 8.314462618
P_ATM = 101325.0

# ==========================================================================
# 1. UNITS
# ==========================================================================
PRESSURE_UNITS = {  # unit: (factor to Pa, gauge?, differential label)
    "barg": (1e5, True, "bar"), "bara": (1e5, False, "bar"),
    "kPag": (1e3, True, "kPa"), "kPaa": (1e3, False, "kPa"),
    "MPag": (1e6, True, "MPa"), "MPaa": (1e6, False, "MPa"),
    "psig": (6894.757, True, "psi"), "psia": (6894.757, False, "psi"),
    "kg/cm2g": (98066.5, True, "kg/cm2"),
}
TEMP_UNITS = ["degC", "degF", "K"]
FLOW_UNITS = {  # unit: (kind, factor to SI)  vol -> m3/s, mass -> kg/s, std -> Nm3/s
    "m3/h": ("vol", 1 / 3600), "L/s": ("vol", 1e-3), "L/min": ("vol", 1e-3 / 60),
    "USgpm": ("vol", 6.30902e-5), "bbl/d": ("vol", 0.158987 / 86400),
    "kg/h": ("mass", 1 / 3600), "kg/s": ("mass", 1.0), "t/h": ("mass", 1000 / 3600),
    "lb/h": ("mass", 0.45359237 / 3600),
    "Nm3/h": ("std0", 1 / 3600), "Sm3/h": ("std15", 1 / 3600),
    "MMSCFD": ("std15", 1e6 * 0.0283168 / 86400 * (273.15 + 15) / (273.15 + 15.56)),
}
LENGTH_UNITS = {"m": 1.0, "km": 1000.0, "ft": 0.3048}


def p_to_pa_abs(val, unit, patm=P_ATM):
    f, gauge, _ = PRESSURE_UNITS[unit]
    return val * f + (patm if gauge else 0.0)


def pa_abs_to_unit(pa, unit, patm=P_ATM):
    f, gauge, _ = PRESSURE_UNITS[unit]
    return (pa - (patm if gauge else 0.0)) / f


def dp_to_unit(dpa, unit):
    return dpa / PRESSURE_UNITS[unit][0]


def dp_label(unit):
    return PRESSURE_UNITS[unit][2]


def t_to_k(val, unit):
    return {"degC": val + 273.15, "degF": (val - 32) * 5 / 9 + 273.15, "K": val}[unit]


def k_to_t(k, unit):
    return {"degC": k - 273.15, "degF": (k - 273.15) * 9 / 5 + 32, "K": k}[unit]


def flow_to_kgs(val, unit, rho_actual, fluid):
    kind, f = FLOW_UNITS[unit]
    if kind == "mass":
        return val * f
    if kind == "vol":
        return val * f * rho_actual
    T_ref = 273.15 if kind == "std0" else 288.15
    return val * f * fluid.props(T_ref, P_ATM)["rho"]


# ==========================================================================
# 2. PIPE DATA  (ASME B36.10M - OD and wall thickness in mm)
# ==========================================================================
PIPE_OD = {"0.5": 21.3, "0.75": 26.7, "1": 33.4, "1.25": 42.2, "1.5": 48.3, "2": 60.3,
           "2.5": 73.0, "3": 88.9, "4": 114.3, "6": 168.3, "8": 219.1, "10": 273.1,
           "12": 323.8, "14": 355.6, "16": 406.4, "18": 457.0, "20": 508.0, "24": 610.0}
WALL = {
    "10S": {"0.5": 2.11, "0.75": 2.11, "1": 2.77, "1.25": 2.77, "1.5": 2.77, "2": 2.77, "2.5": 3.05,
            "3": 3.05, "4": 3.05, "6": 3.40, "8": 3.76, "10": 4.19, "12": 4.57, "14": 4.78,
            "16": 4.78, "18": 4.78, "20": 5.54, "24": 6.35},
    "40":  {"0.5": 2.77, "0.75": 2.87, "1": 3.38, "1.25": 3.56, "1.5": 3.68, "2": 3.91, "2.5": 5.16,
            "3": 5.49, "4": 6.02, "6": 7.11, "8": 8.18, "10": 9.27, "12": 10.31, "14": 11.13,
            "16": 12.70, "18": 14.27, "20": 15.09, "24": 17.48},
    "STD": {"0.5": 2.77, "0.75": 2.87, "1": 3.38, "1.25": 3.56, "1.5": 3.68, "2": 3.91, "2.5": 5.16,
            "3": 5.49, "4": 6.02, "6": 7.11, "8": 8.18, "10": 9.27, "12": 9.53, "14": 9.53,
            "16": 9.53, "18": 9.53, "20": 9.53, "24": 9.53},
    "80":  {"0.5": 3.73, "0.75": 3.91, "1": 4.55, "1.25": 4.85, "1.5": 5.08, "2": 5.54, "2.5": 7.01,
            "3": 7.62, "4": 8.56, "6": 10.97, "8": 12.70, "10": 15.09, "12": 17.48, "14": 19.05,
            "16": 21.44, "18": 23.83, "20": 26.19, "24": 30.96},
    "XS":  {"0.5": 3.73, "0.75": 3.91, "1": 4.55, "1.25": 4.85, "1.5": 5.08, "2": 5.54, "2.5": 7.01,
            "3": 7.62, "4": 8.56, "6": 10.97, "8": 12.70, "10": 12.70, "12": 12.70, "14": 12.70,
            "16": 12.70, "18": 12.70, "20": 12.70, "24": 12.70},
}
MATERIALS = {  # roughness [mm], wall conductivity [W/mK]
    "Carbon steel (new)":      {"rough_mm": 0.045, "k": 45.0},
    "Carbon steel (corroded)": {"rough_mm": 0.15, "k": 45.0},
    "Stainless steel":         {"rough_mm": 0.015, "k": 16.0},
    "Galvanised steel":        {"rough_mm": 0.15, "k": 45.0},
    "Cast / ductile iron":     {"rough_mm": 0.26, "k": 50.0},
    "Copper":                  {"rough_mm": 0.0015, "k": 380.0},
    "HDPE":                    {"rough_mm": 0.007, "k": 0.45},
    "PVC":                     {"rough_mm": 0.0015, "k": 0.19},
}


def size_dims(size):
    """size = {'nps': '6', 'sch': '40'} or {'nps': 'custom', 'id_mm': 150, 'wall_mm': 6}.
    Returns (ID [m], OD [m], label)."""
    nps = str(size.get("nps", "custom"))
    if nps == "custom":
        idm = float(size.get("id_mm", 0) or 0)
        wall = float(size.get("wall_mm", 6) or 6)
        if idm <= 0:
            raise ValueError("Custom pipe ID must be > 0 mm")
        return idm / 1000, (idm + 2 * wall) / 1000, f"ID {idm:g} mm"
    sch = str(size.get("sch", "40"))
    od = PIPE_OD[nps]
    wt = WALL[sch][nps]
    return (od - 2 * wt) / 1000, od / 1000, f'{nps}" Sch {sch}'


def size_key(size):
    """Used for mismatch detection (nominal size only)."""
    nps = str(size.get("nps", "custom"))
    if nps == "custom":
        return f"ID{round(float(size.get('id_mm', 0) or 0), 1)}"
    return f"NPS{nps}"


# ==========================================================================
# 3. FLUIDS
# ==========================================================================
def _sutherland(mu0, T0, S):
    return lambda T: mu0 * (T / T0) ** 1.5 * (T0 + S) / (T + S)


def _water_rho(T):
    c = T - 273.15
    return 1000 * (1 - (c + 288.9414) / (508929.2 * (c + 68.12963)) * (c - 3.9863) ** 2)


def _water_mu(T):
    return 2.414e-5 * 10 ** (247.8 / (T - 140))


def _water_k(T):
    r = T / 298.15
    return 0.6065 * (-1.48445 + 4.12292 * r - 1.63866 * r * r)


def water_psat(T):
    """Antoine, Pa (valid ~1-370 C)."""
    c = T - 273.15
    if c <= 100:
        a, b, cc = 8.07131, 1730.63, 233.426
    else:
        a, b, cc = 8.14019, 1810.94, 244.485
    return 10 ** (a - b / (cc + c)) * 133.322


def water_tsat(P):
    lo, hi = 274.0, 640.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if water_psat(mid) > P:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _andrade(rho15, nu1_cst, T1c, nu2_cst, T2c):
    """Fit ln(mu) = A + B/T through two kinematic viscosity points."""
    T1, T2 = T1c + 273.15, T2c + 273.15
    rho1 = rho15 * (1 - 0.00072 * (T1c - 15))
    rho2 = rho15 * (1 - 0.00072 * (T2c - 15))
    m1, m2 = nu1_cst * 1e-6 * rho1, nu2_cst * 1e-6 * rho2
    B = math.log(m1 / m2) / (1 / T1 - 1 / T2)
    A = math.log(m1) - B / T1
    return lambda T: math.exp(A + B / T)


class Fluid:
    def __init__(self, name, phase, rho=None, mu=None, cp=None, k=None, mw=None, gamma=1.4,
                 z=1.0, mu_jt=0.0, psat=None, note=""):
        self.name, self.phase, self.note = name, phase, note
        self._rho, self._mu, self._cp, self._k = rho, mu, cp, k
        self.mw, self.gamma, self.z, self.mu_jt, self._psat = mw, gamma, z, mu_jt, psat

    def props(self, T, P):
        if self.phase == "gas":
            rho = P * self.mw / 1000 / (self.z * R_UNIV * T)
            c = math.sqrt(self.gamma * self.z * R_UNIV * T / (self.mw / 1000))
        else:
            rho = self._rho(T) if callable(self._rho) else self._rho
            c = 1400.0
        val = lambda f: f(T) if callable(f) else f
        return {"rho": rho, "mu": val(self._mu), "cp": val(self._cp), "k": val(self._k),
                "c": c, "mu_jt": self.mu_jt,
                "psat": self._psat(T) if self._psat else None}

    def check_phase(self, T, P):
        """Return a warning keyword if the built-in single-phase model is violated."""
        if self.name.startswith("Water") and self.phase == "liquid" and P < water_psat(T):
            return "flash"
        if self.name.startswith("Steam") and T < water_tsat(P):
            return "condense"
        return None


def builtin_fluids():
    F = {}
    F["Water"] = Fluid("Water", "liquid", _water_rho, _water_mu, 4184.0, _water_k, psat=water_psat)
    F["Sea water (approx.)"] = Fluid("Sea water (approx.)", "liquid", lambda T: 1.025 * _water_rho(T),
                                     lambda T: 1.08 * _water_mu(T), 3990.0, lambda T: 0.98 * _water_k(T),
                                     psat=water_psat)
    F["MEG 30% (approx.)"] = Fluid("MEG 30% (approx.)", "liquid", lambda T: 1045 - 0.45 * (T - 293.15),
                                   _andrade(1048, 2.4, 20, 0.95, 80), 3700.0, 0.48)
    F["Diesel"] = Fluid("Diesel", "liquid", lambda T: 840 * (1 - 0.00072 * (T - 288.15)),
                        _andrade(840, 3.0, 40, 1.2, 100), 2000.0, 0.13)
    F["Light crude oil"] = Fluid("Light crude oil", "liquid", lambda T: 850 * (1 - 0.00070 * (T - 288.15)),
                                 _andrade(850, 6.0, 40, 2.2, 100), 1950.0, 0.13)
    F["Heavy fuel oil"] = Fluid("Heavy fuel oil", "liquid", lambda T: 980 * (1 - 0.00064 * (T - 288.15)),
                                _andrade(980, 380.0, 50, 30.0, 100), 1800.0, 0.12)
    F["Air"] = Fluid("Air", "gas", mu=_sutherland(1.716e-5, 273.15, 110.4), cp=1005.0,
                     k=lambda T: 0.0241 * (T / 273.15) ** 0.81, mw=28.96, gamma=1.40, mu_jt=0.25e-5)
    F["Nitrogen"] = Fluid("Nitrogen", "gas", mu=_sutherland(1.663e-5, 273.15, 107.0), cp=1040.0,
                          k=lambda T: 0.0240 * (T / 273.15) ** 0.80, mw=28.013, gamma=1.40, mu_jt=0.22e-5)
    F["Natural gas (CH4-rich)"] = Fluid("Natural gas (CH4-rich)", "gas", mu=_sutherland(1.03e-5, 273.15, 198.0),
                                         cp=2220.0, k=lambda T: 0.0302 * (T / 273.15) ** 1.2, mw=17.4,
                                         gamma=1.30, mu_jt=0.45e-5)
    F["Carbon dioxide"] = Fluid("Carbon dioxide", "gas", mu=_sutherland(1.37e-5, 273.15, 222.0), cp=845.0,
                                k=lambda T: 0.0146 * (T / 273.15) ** 1.3, mw=44.01, gamma=1.29, mu_jt=1.1e-5)
    F["Hydrogen"] = Fluid("Hydrogen", "gas", mu=_sutherland(8.41e-6, 273.15, 72.0), cp=14300.0,
                          k=lambda T: 0.168 * (T / 273.15) ** 0.8, mw=2.016, gamma=1.41, mu_jt=-0.03e-5)
    F["Steam (superheated)"] = Fluid("Steam (superheated)", "gas",
                                     mu=lambda T: 1.23e-5 + 4.0e-8 * (T - 373.15), cp=2080.0,
                                     k=lambda T: 0.025 + 8.0e-5 * (T - 373.15), mw=18.015, gamma=1.30,
                                     z=0.97, mu_jt=1.5e-5)
    return F


# ---- Optional CoolProp (uncomment in requirements.txt for real-fluid properties)
try:
    from CoolProp.CoolProp import PropsSI, PhaseSI
    HAS_COOLPROP = True
except Exception:  # pragma: no cover
    HAS_COOLPROP = False

COOLPROP_FLUIDS = ["Water", "Nitrogen", "Air", "Methane", "CarbonDioxide", "Ammonia",
                   "Hydrogen", "Oxygen", "Propane", "R134a"]


def _cp_phase(T, P, name):
    """'liquid' or 'gas' from CoolProp PhaseSI ('twophase' -> None)."""
    ph = PhaseSI("T", T, "P", P, name).lower()
    if "twophase" in ph:
        return None
    return "liquid" if "liquid" in ph else "gas"


class CPFluid(Fluid):
    """Real-fluid properties via CoolProp (optional dependency)."""
    def __init__(self, cp_name):
        super().__init__("CoolProp: " + cp_name, "liquid")
        self.cp_name = cp_name

    def set_phase_from(self, T, P):
        ph = _cp_phase(T, P, self.cp_name)
        if ph is None:
            raise ValueError("Inlet state is two-phase - this tool handles single-phase flow only.")
        self.phase = ph

    def props(self, T, P):
        n = self.cp_name
        d = {"rho": PropsSI("D", "T", T, "P", P, n), "mu": PropsSI("V", "T", T, "P", P, n),
             "cp": PropsSI("C", "T", T, "P", P, n), "k": PropsSI("L", "T", T, "P", P, n)}
        try:
            d["c"] = PropsSI("A", "T", T, "P", P, n)
        except Exception:
            d["c"] = 1400.0
        try:
            d["mu_jt"] = PropsSI("d(T)/d(P)|Hmass", "T", T, "P", P, n)
        except Exception:
            d["mu_jt"] = 0.0
        d["psat"] = None
        if self.phase == "liquid":
            try:
                d["psat"] = PropsSI("P", "T", T, "Q", 0, n)
            except Exception:
                pass
        return d

    def check_phase(self, T, P):
        try:
            ph = _cp_phase(T, P, self.cp_name)
        except Exception:
            ph = None
        if ph != self.phase:
            return "flash" if self.phase == "liquid" else "condense"
        return None


def make_fluid(inlet):
    name = inlet.get("fluid", "Water")
    if name.startswith("CoolProp: ") and HAS_COOLPROP:
        f = CPFluid(name.replace("CoolProp: ", ""))
        f.set_phase_from(t_to_k(float(inlet["T"]), inlet["T_unit"]),
                         p_to_pa_abs(float(inlet["P"]), inlet["P_unit"]))
        return f
    if name == "Custom liquid":
        c = inlet.get("custom", {})
        return Fluid(name, "liquid", float(c.get("rho", 1000)), float(c.get("mu_cP", 1)) / 1000,
                     float(c.get("cp", 4184)), float(c.get("k", 0.6)))
    if name == "Custom gas":
        c = inlet.get("custom", {})
        return Fluid(name, "gas", mu=float(c.get("mu_cP", 0.018)) / 1000, cp=float(c.get("cp", 1005)),
                     k=float(c.get("k", 0.026)), mw=float(c.get("mw", 28.96)),
                     gamma=float(c.get("gamma", 1.4)), z=float(c.get("z", 1.0)),
                     mu_jt=float(c.get("mu_jt_K_per_bar", 0.25)) * 1e-5)
    F = builtin_fluids()
    if name not in F:
        raise ValueError(f"Unknown fluid '{name}'")
    return F[name]


# ==========================================================================
# 4. FRICTION & FITTINGS
# ==========================================================================
def colebrook(Re, rel_rough):
    """Darcy friction factor, Colebrook-White (iterated from Swamee-Jain seed)."""
    f = 0.25 / math.log10(rel_rough / 3.7 + 5.74 / Re ** 0.9) ** 2
    for _ in range(30):
        f_new = (-2.0 * math.log10(rel_rough / 3.7 + 2.51 / (Re * math.sqrt(f)))) ** -2
        if abs(f_new - f) < 1e-10:
            return f_new
        f = f_new
    return f


def darcy_f(Re, rel_rough):
    if Re <= 0:
        return 0.0
    if Re < 2300:
        return 64.0 / Re
    if Re < 4000:  # transitional - linear blend
        w = (Re - 2300) / 1700
        return (1 - w) * 64 / 2300 + w * colebrook(4000, rel_rough)
    return colebrook(Re, rel_rough)


def f_turbulent(D, rough):
    """Fully turbulent friction factor f_T used by Crane K = n*f_T."""
    return 0.25 / math.log10(rough / (3.7 * D)) ** 2


def regime(Re):
    return "Laminar" if Re < 2300 else ("Transitional" if Re < 4000 else "Turbulent")


# Crane TP-410 style: K = n * f_T
ELBOWS = {"lr90": ("90 deg LR elbow (r/D=1.5)", 14), "sr90": ("90 deg SR elbow (r/D=1)", 30),
          "mitre90": ("90 deg mitre elbow", 60)}
FITTINGS = {
    "gate":      ("Gate valve (open)", 8),
    "ball":      ("Ball valve (full bore)", 3),
    "globe":     ("Globe valve (open)", 340),
    "butterfly": ("Butterfly valve (open)", 45),
    "plug":      ("Plug valve (straight)", 18),
    "check":     ("Swing check valve", 100),
    "liftcheck": ("Lift check valve", 600),
    "elbow45":   ("45 deg elbow (in-line)", 16),
    "strainer":  ("Y-strainer (clean, approx.)", 150),
    "customK":   ("Custom K value", None),
    "fixeddp":   ("Fixed dP (equipment / meter)", None),
}
TEE_K = {"run": 20, "branch": 60}


def reducer_K(D_in, D_out, angle_deg):
    """Crane formulae, K referenced to the SMALL-diameter velocity."""
    small, large = min(D_in, D_out), max(D_in, D_out)
    beta = small / large
    th = math.radians(max(1.0, min(180.0, angle_deg)))
    if D_out < D_in:  # contraction
        K = 0.8 * math.sin(th / 2) * (1 - beta ** 2) if angle_deg <= 45 \
            else 0.5 * math.sqrt(math.sin(th / 2)) * (1 - beta ** 2)
    else:  # expansion
        K = 2.6 * math.sin(th / 2) * (1 - beta ** 2) ** 2 if angle_deg <= 45 else (1 - beta ** 2) ** 2
    return K, small


# ==========================================================================
# 5. 3D GEOMETRY (same rules as index.html)
# ==========================================================================
# E = +x, N = +y, Up = +z
DIRS = {"east": (1, 0, 0), "west": (-1, 0, 0), "north": (0, 1, 0), "south": (0, -1, 0),
        "up": (0, 0, 1), "down": (0, 0, -1)}
LEGACY_DIR = {"right": "east", "left": "west"}  # files saved by the old 2D version
PERP = {d: [k for k, v in DIRS.items() if sum(a * b for a, b in zip(v, DIRS[d])) == 0] for d in DIRS}


def normalize(net):
    """Map legacy 2D directions to 3D compass directions (in place)."""
    inl = net.get("inlet", {})
    hd = inl.get("heading", "east")
    inl["heading"] = LEGACY_DIR.get(hd, hd)
    for s in net.get("segments", []):
        for k in ("turn", "branch_dir"):
            if k in s:
                s[k] = LEGACY_DIR.get(s[k], s[k])
    return net


def out_ports(seg_type):
    return {"tee": ["run", "branch"], "outlet": []}.get(seg_type, ["out"])


# ==========================================================================
# 6. VALIDATION
# ==========================================================================
def validate(net):
    errors = []
    segs = net.get("segments", [])
    if not segs:
        return ["Add at least one segment after the inlet."]
    by_id = {s["id"]: s for s in segs}
    roots = [s for s in segs if s.get("parent") in (None, "", "INLET")]
    if len(roots) != 1:
        errors.append("Exactly one segment must be connected to the inlet.")
    used = set()
    for s in segs:
        p = s.get("parent")
        if p in (None, "", "INLET"):
            key = ("INLET", "out")
        else:
            if p not in by_id:
                errors.append(f"{s['id']}: parent {p} not found.")
                continue
            if s.get("port") not in out_ports(by_id[p]["type"]):
                errors.append(f"{s['id']}: invalid port '{s.get('port')}' on {p}.")
            key = (p, s.get("port"))
        if key in used:
            errors.append(f"{s['id']}: port {key[1]} of {key[0]} already has a segment.")
        used.add(key)
    for s in segs:
        for port in out_ports(s["type"]):
            if (s["id"], port) not in used:
                errors.append(f"Branch end at {s['id']} ({port}) is open - terminate it with an Outlet.")
    children = {}
    for s in segs:
        children.setdefault(s.get("parent") or "INLET", []).append(s)
    try:
        inlet_size = net["inlet"]["size"]
    except KeyError:
        return errors + ["Inlet pipe size missing."]

    def walk(parent_id, up_size, heading):
        for c in children.get(parent_id, []):
            t = c["type"]
            port = c.get("port", "out")
            h = heading
            if parent_id != "INLET" and by_id[parent_id]["type"] == "tee" and port == "branch":
                h = by_id[parent_id].get("branch_dir")
            if t == "outlet":
                continue
            this_size = c.get("size", up_size)
            if t == "reducer":
                if size_key(this_size) == size_key(up_size):
                    errors.append(f"{c['id']}: reducer outlet size equals inlet size.")
            elif size_key(this_size) != size_key(up_size):
                errors.append(f"SIZE MISMATCH at {c['id']}: {size_dims(this_size)[2]} connected to "
                              f"{size_dims(up_size)[2]}. Insert a reducer/expander.")
            if t == "elbow" and c.get("turn") not in PERP.get(h, []):
                errors.append(f"{c['id']}: elbow cannot turn '{c.get('turn')}' when flow is heading '{h}'.")
            if t == "tee" and c.get("branch_dir") not in PERP.get(h, []):
                errors.append(f"{c['id']}: tee branch '{c.get('branch_dir')}' must be perpendicular to '{h}'.")
            new_h = c.get("turn") if t == "elbow" else h
            walk(c["id"], this_size, new_h)

    walk("INLET", inlet_size, net["inlet"].get("heading", "east"))
    return errors


# ==========================================================================
# 7. SOLVER
# ==========================================================================
def solve(net):
    out = {"ok": False, "errors": [], "warnings": [], "infeasible": [], "rows": [], "profiles": [],
           "nodes": {}, "outlets": [], "summary": {}}
    net = normalize(copy.deepcopy(net))
    try:
        errs = validate(net)
        if errs:
            out["errors"] = errs
            return out
        return _solve(net, out)
    except Exception as e:  # report instead of crashing Streamlit
        out["errors"].append(f"Calculation error: {e}")
        return out


def _solve(net, out):
    inlet = net["inlet"]
    opt = net.get("options", {})
    patm = float(opt.get("P_atm_bar", 1.01325)) * 1e5
    Pu, Tu, Qu = inlet["P_unit"], inlet["T_unit"], inlet["Q_unit"]
    dPu = dp_label(Pu)
    lu = opt.get("len_unit", "m")
    Lmul = LENGTH_UNITS.get(lu, 1.0)

    fluid = make_fluid(inlet)
    P0 = p_to_pa_abs(float(inlet["P"]), Pu, patm)
    T0 = t_to_k(float(inlet["T"]), Tu)
    if P0 <= 0:
        out["errors"].append("Inlet absolute pressure must be > 0.")
        return out
    pr0 = fluid.props(T0, P0)
    m_in = flow_to_kgs(float(inlet["Q"]), Qu, pr0["rho"], fluid)
    if m_in <= 0:
        out["errors"].append("Inlet flow must be > 0.")
        return out

    mat = MATERIALS.get(net.get("pipe", {}).get("material", "Carbon steel (new)"), MATERIALS["Carbon steel (new)"])
    rough = float(net.get("pipe", {}).get("rough_mm") or mat["rough_mm"]) / 1000
    k_wall = mat["k"]

    th = net.get("thermal", {})
    Tamb = t_to_k(float(th.get("T_amb", 30)), Tu)
    tcfg = thermal.thermal_config(th, Tamb)

    segs = net["segments"]
    by_id = {s["id"]: s for s in segs}
    children = {}
    for s in segs:
        children.setdefault((s.get("parent") or "INLET", s.get("port", "out")), []).append(s)

    # ---- outlet flows (mass balance)
    outlets = [s for s in segs if s["type"] == "outlet"]
    blanks = [o for o in outlets if o.get("Q") in (None, "")]
    if len(outlets) == 1 and blanks:
        outlets[0]["_m"] = m_in
    else:
        spec = 0.0
        for o in outlets:
            if o not in blanks:
                o["_m"] = flow_to_kgs(float(o["Q"]), o.get("Q_unit", Qu), pr0["rho"], fluid)
                spec += o["_m"]
        if len(blanks) > 1:
            out["errors"].append("Only one outlet flow may be left blank (it takes the remainder).")
            return out
        if blanks:
            rem = m_in - spec
            if rem <= 0:
                out["errors"].append(f"Outlet flows ({spec*3600:.1f} kg/h) exceed inlet flow ({m_in*3600:.1f} kg/h).")
                return out
            blanks[0]["_m"] = rem
        elif abs(spec - m_in) / m_in > 0.005:
            out["errors"].append(f"Mass balance error: sum of outlet flows = {spec*3600:.1f} kg/h, "
                                 f"inlet = {m_in*3600:.1f} kg/h. Adjust outlet flows (or leave one blank).")
            return out
        for o in outlets:
            if o["_m"] <= 0:
                out["errors"].append(f"Outlet {o.get('label', o['id'])}: flow must be > 0.")
                return out

    mflow = {}

    def subtree_m(sid):
        s = by_id[sid]
        if s["type"] == "outlet":
            mflow[sid] = s["_m"]
            return s["_m"]
        tot = 0.0
        for port in out_ports(s["type"]):
            for c in children.get((sid, port), []):
                tot += subtree_m(c["id"])
        mflow[sid] = tot
        return tot

    root = children[("INLET", "out")][0]
    subtree_m(root["id"])

    def port_m(sid, port):
        return sum(mflow[c["id"]] for c in children.get((sid, port), []))

    def disp_P(p):
        return pa_abs_to_unit(p, Pu, patm)

    def disp_T(t):
        return k_to_t(t, Tu)

    flags = {"infeasible": False}

    def check_state(sid, label, P, T, v, pr):
        if P <= 1000:
            out["infeasible"].append(f"dP NOT FEASIBLE: pressure at {label} falls to "
                                     f"{disp_P(max(P, 0)):.3f} {Pu} (~vacuum). Available pressure is "
                                     f"insufficient for this flow / routing.")
            flags["infeasible"] = True
            return False
        if fluid.phase == "liquid" and pr.get("psat") and P < pr["psat"]:
            out["infeasible"].append(f"dP NOT FEASIBLE: pressure at {label} ({disp_P(P):.3f} {Pu}) is below "
                                     f"vapour pressure ({disp_P(pr['psat']):.3f} {Pu}) -> flashing/cavitation.")
            flags["infeasible"] = True
            return False
        if fluid.phase == "gas" and pr["c"] > 0:
            mach = v / pr["c"]
            if mach >= 1.0:
                out["infeasible"].append(f"dP NOT FEASIBLE: choked flow at {label} (Mach {mach:.2f}).")
                flags["infeasible"] = True
                return False
            if mach > 0.5:
                out["warnings"].append(f"{sid}: high Mach number {mach:.2f} (> 0.5) - compressible effects significant.")
        if fluid.check_phase(T, P) == "condense":
            out["warnings"].append(f"{sid}: temperature below saturation at {label} - condensation likely.")
        return True

    vel_warned = set()

    def check_velocity(sid, v):
        if sid in vel_warned:
            return
        if fluid.phase == "liquid":
            if v > 3.0:
                out["warnings"].append(f"{sid}: liquid velocity {v:.2f} m/s > 3 m/s (erosion / noise / surge risk).")
                vel_warned.add(sid)
            elif v < 0.5:
                out["warnings"].append(f"{sid}: liquid velocity {v:.2f} m/s < 0.5 m/s (possible settling / oversized).")
                vel_warned.add(sid)
        else:
            lim = 40.0 if "Steam" in fluid.name else 30.0
            if v > lim:
                out["warnings"].append(f"{sid}: gas velocity {v:.1f} m/s > {lim:.0f} m/s (noise / erosion).")
                vel_warned.add(sid)

    rows = out["rows"]

    def add_row(seg, kind, label, D, L, dz, m, st_in, st_out, v1, v2, Re, f, K, dPf, dPk, dPz, dPa, Q_loss, size_lbl):
        pr = st_in["pr"]
        rows.append({
            "Segment": seg["id"], "Type": kind, "Description": label, "Size": size_lbl,
            "ID [mm]": round(D * 1000, 1), f"L [{lu}]": round(L / Lmul, 3),
            "dz [m]": round(dz, 3),
            "Mass flow [kg/h]": round(m * 3600, 2), "Q actual [m3/h]": round(m / pr["rho"] * 3600, 3),
            "Density [kg/m3]": round(pr["rho"], 3), "Viscosity [cP]": round(pr["mu"] * 1000, 4),
            "Velocity in [m/s]": round(v1, 3), "Velocity out [m/s]": round(v2, 3),
            "Re [-]": round(Re), "Regime": regime(Re),
            "Darcy f [-]": round(f, 5) if f else None, "K [-]": round(K, 3) if K else None,
            f"dP friction [{dPu}]": dp_to_unit(dPf, Pu), f"dP minor [{dPu}]": dp_to_unit(dPk, Pu),
            f"dP static [{dPu}]": dp_to_unit(dPz, Pu), f"dP kinetic [{dPu}]": dp_to_unit(dPa, Pu),
            f"dP total [{dPu}]": dp_to_unit(st_in["P"] - st_out["P"], Pu),
            f"P in [{Pu}]": disp_P(st_in["P"]), f"P out [{Pu}]": disp_P(st_out["P"]),
            f"T in [{Tu}]": disp_T(st_in["T"]), f"T out [{Tu}]": disp_T(st_out["T"]),
            "Heat loss [kW]": round(Q_loss / 1000, 3),
        })

    def set_xyz(row, pt):
        row["End E [m]"] = round(pt["E"], 2)
        row["End N [m]"] = round(pt["N"], 2)
        row["End Z [m]"] = round(pt["z"], 2)

    def pipe(seg, st, D, Do, m, heading, size_lbl, path):
        L = float(seg.get("L", 0) or 0) * Lmul
        if L <= 0:
            raise ValueError(f"{seg['id']}: pipe length must be > 0")
        vec = DIRS[heading]
        vertical = heading in ("up", "down")
        if vertical:
            dz_tot = vec[2] * L
        else:
            dz_tot = float(seg.get("dz", 0) or 0)
            if abs(dz_tot) >= L:
                raise ValueError(f"{seg['id']}: elevation change |dz| must be smaller than the pipe length")
        Lh = math.sqrt(max(L * L - dz_tot * dz_tot, 0.0))  # horizontal projection
        A = math.pi * D * D / 4
        n = int(min(400 if fluid.phase == "gas" else 200, max(10, math.ceil(L / 10))))
        dL, dzs = L / n, dz_tot / n
        dE, dN = vec[0] * Lh / n, vec[1] * Lh / n
        P, T = st["P"], st["T"]
        sums = dict(f=0.0, z=0.0, a=0.0, q=0.0)
        Re_list, f_list = [], []
        pr_in = fluid.props(T, P)
        v_first = m / (pr_in["rho"] * A)
        v = v_first
        for i in range(n):
            pr1 = fluid.props(T, P)
            v1 = m / (pr1["rho"] * A)                                            # continuity
            Re = pr1["rho"] * v1 * D / pr1["mu"]                                 # Reynolds
            f = darcy_f(Re, rough / D)                                           # friction factor
            P2, rho_avg = P, pr1["rho"]
            for _ in range(3 if fluid.phase == "gas" else 1):
                dPf = f * dL / D * rho_avg * (m / (rho_avg * A)) ** 2 / 2       # Darcy-Weisbach
                dPz = rho_avg * G * dzs                                          # static head
                rho2 = fluid.props(T, max(P2, 1000))["rho"] if fluid.phase == "gas" else pr1["rho"]
                dPa = (m / A) ** 2 * (1 / rho2 - 1 / pr1["rho"])                 # acceleration
                P2 = P - dPf - dPz - dPa                                         # energy balance
                rho_avg = 0.5 * (pr1["rho"] + rho2)
            T2, q = T, 0.0
            if tcfg["enabled"]:
                Pr_num = pr1["cp"] * pr1["mu"] / pr1["k"]
                hi = thermal.inside_htc(Re, Pr_num, pr1["k"], D, f)
                R = thermal.resistance_per_length(T, hi, D, Do, k_wall, tcfg)["R"]
                T2, q = thermal.heat_loss_step(T, m, pr1["cp"], R, dL, Tamb)
            T2 += thermal.pressure_effect_dT(fluid.phase, P2 - P, dPf, pr1["rho"], pr1["cp"], pr1["mu_jt"])
            if P2 <= 1000:  # pressure exhausted - clamp for reporting, check_state flags it
                P2 = 0.0
            sums["f"] += dPf; sums["z"] += dPz; sums["a"] += dPa; sums["q"] += q
            Re_list.append(Re); f_list.append(f)
            P, T = P2, T2
            v = m / (fluid.props(T, max(P, 1000))["rho"] * A) if P > 1000 else v1
            last = path[-1]
            path.append({"x": last["x"] + dL, "P": P, "T": T, "E": last["E"] + dE, "N": last["N"] + dN,
                         "z": last["z"] + dzs, "seg": seg["id"]})
            if not check_state(seg["id"], f"{seg['id']} (+{(i + 1) * dL / Lmul:.1f} {lu})",
                               P, T, v, fluid.props(T, max(P, 1000))):
                break
        check_velocity(seg["id"], v_first)
        Re_m, f_m = sum(Re_list) / len(Re_list), sum(f_list) / len(f_list)
        st2 = {"P": P, "T": T, "pr": fluid.props(T, max(P, 1000))}
        if heading == "up":
            kind = "Pipe - Vertical up"
        elif heading == "down":
            kind = "Pipe - Vertical down"
        elif dz_tot == 0:
            kind = "Pipe - Horizontal"
        else:
            kind = "Pipe - Sloped " + ("rising" if dz_tot > 0 else "falling")
        desc = f"{L / Lmul:g} {lu} {heading}" + (f", dz {dz_tot:+.2f} m" if (not vertical and dz_tot) else "")
        add_row(seg, kind, desc, D, L, dz_tot, m, {"P": st["P"], "T": st["T"], "pr": pr_in}, st2,
                v_first, v, Re_m, f_m, None, sums["f"], 0.0, sums["z"], sums["a"], sums["q"], size_lbl)
        set_xyz(rows[-1], path[-1])
        return st2

    def local_loss(seg, st, D_in, D_out, K, m, kind, label, size_lbl, path, fixed_dp=0.0):
        """Zero-length element: P_out = P_in - K rho v^2/2 - (rho/2)(v_out^2 - v_in^2) - fixed."""
        pr = fluid.props(st["T"], st["P"])
        A_in, A_out = math.pi * D_in ** 2 / 4, math.pi * D_out ** 2 / 4
        v_in, v_out = m / (pr["rho"] * A_in), m / (pr["rho"] * A_out)
        D_ref = min(D_in, D_out)
        v_ref = m / (pr["rho"] * math.pi * D_ref ** 2 / 4)
        Re = pr["rho"] * v_ref * D_ref / pr["mu"]
        dPk = K * pr["rho"] * v_ref ** 2 / 2 + fixed_dp
        dPa = pr["rho"] * (v_out ** 2 - v_in ** 2) / 2
        P2 = st["P"] - dPk - dPa
        T2 = st["T"] + thermal.pressure_effect_dT(fluid.phase, P2 - st["P"], dPk, pr["rho"], pr["cp"], pr["mu_jt"])
        st2 = {"P": P2, "T": T2, "pr": fluid.props(T2, max(P2, 1000))}
        add_row(seg, kind, label, D_ref, 0.0, 0.0, m, {"P": st["P"], "T": st["T"], "pr": pr}, st2,
                v_in, v_out, Re, f_turbulent(D_ref, rough), K, 0.0, dPk, 0.0, dPa, 0.0, size_lbl)
        path.append(dict(path[-1], P=P2, T=T2, seg=seg["id"]))
        set_xyz(rows[-1], path[-1])
        check_velocity(seg["id"], v_ref)
        check_state(seg["id"], seg["id"], P2, T2, v_ref, st2["pr"])
        return st2

    def visit(seg, st, size, heading, path):
        if flags["infeasible"]:
            return
        sid, t = seg["id"], seg["type"]
        D, Do, lbl = size_dims(seg.get("size", size)) if t != "outlet" else size_dims(size)
        m = mflow[sid]
        out["nodes"][sid] = {"P_in": disp_P(st["P"]), "T_in": disp_T(st["T"])}

        if t == "outlet":
            Pmin = seg.get("Pmin")
            P_disp = disp_P(st["P"])
            rec = {"id": sid, "label": seg.get("label") or sid, "P": P_disp, "T": disp_T(st["T"]),
                   "m_kgh": m * 3600, "Q_m3h": m / st["pr"]["rho"] * 3600,
                   "dP_total": dp_to_unit(P0 - st["P"], Pu), "length": path[-1]["x"] / Lmul,
                   "Pmin": Pmin, "feasible": True,
                   "E": path[-1]["E"], "N": path[-1]["N"], "Z": path[-1]["z"]}
            if Pmin not in (None, ""):
                req = p_to_pa_abs(float(Pmin), Pu, patm)
                if st["P"] < req:
                    rec["feasible"] = False
                    flags["infeasible"] = True
                    out["infeasible"].append(
                        f"dP NOT FEASIBLE at outlet '{rec['label']}': arrives at {P_disp:.3f} {Pu} but needs "
                        f">= {float(Pmin):.3f} {Pu}. Available dP = {dp_to_unit(P0 - req, Pu):.3f} {dPu}, "
                        f"required dP = {dp_to_unit(P0 - st['P'], Pu):.3f} {dPu}.")
            out["outlets"].append(rec)
            out["profiles"].append({"outlet": rec["label"], "points": [
                {"x": p["x"] / Lmul, "P": disp_P(p["P"]), "T": disp_T(p["T"]), "z": p["z"],
                 "E": p["E"], "N": p["N"], "seg": p["seg"]} for p in path]})
            out["nodes"][sid]["P_out"] = P_disp
            return

        if t == "pipe":
            st2 = pipe(seg, st, D, Do, m, heading, lbl, path)
            new_h = heading
        elif t == "elbow":
            fT = f_turbulent(D, rough)
            name, n = ELBOWS.get(seg.get("elbow", "lr90"), ELBOWS["lr90"])
            st2 = local_loss(seg, st, D, D, n * fT, m, f"Elbow -> {seg.get('turn')}", f"{name}, K={n}fT",
                             lbl, path)
            new_h = seg.get("turn")
        elif t == "fitting":
            fT = f_turbulent(D, rough)
            ftype = seg.get("fitting", "gate")
            name, n = FITTINGS.get(ftype, FITTINGS["gate"])
            fixed = 0.0
            if ftype == "customK":
                K = float(seg.get("K", 0) or 0)
            elif ftype == "fixeddp":
                K = 0.0
                fixed = float(seg.get("dp", 0) or 0) * PRESSURE_UNITS[Pu][0]
                name = f"Fixed dP {seg.get('dp')} {dPu}"
            else:
                K = n * fT
            st2 = local_loss(seg, st, D, D, K, m, "Fitting / valve", name, lbl, path, fixed)
            new_h = heading
        elif t == "reducer":
            D_up = size_dims(size)[0]
            K, _ = reducer_K(D_up, D, float(seg.get("angle", 30) or 30))
            kind = "Reducer" if D < D_up else "Expander"
            st2 = local_loss(seg, st, D_up, D, K, m, kind, f"{size_dims(size)[2]} -> {lbl}", lbl, path)
            new_h = heading
        elif t == "tee":
            fT = f_turbulent(D, rough)
            out["nodes"][sid]["ports"] = {}
            for port in ("run", "branch"):
                mp = port_m(sid, port)
                if mp <= 0:
                    continue
                p_path = list(path)
                st_p = local_loss({"id": f"{sid}.{port}"}, st, D, D, TEE_K[port] * fT, mp,
                                  f"Tee - {port}", f"Equal tee, flow through {port}, K={TEE_K[port]}fT",
                                  lbl, p_path)
                # Bernoulli recovery: tee inlet velocity is based on the total flow
                pr = fluid.props(st["T"], st["P"])
                A = math.pi * D * D / 4
                v_tot, v_p = m / (pr["rho"] * A), mp / (pr["rho"] * A)
                rec = pr["rho"] * (v_tot ** 2 - v_p ** 2) / 2
                st_p["P"] += rec
                p_path[-1]["P"] = st_p["P"]
                rows[-1][f"dP kinetic [{dPu}]"] = dp_to_unit(-rec, Pu)
                rows[-1][f"dP total [{dPu}]"] = dp_to_unit(st["P"] - st_p["P"], Pu)
                rows[-1][f"P out [{Pu}]"] = disp_P(st_p["P"])
                rows[-1]["Velocity in [m/s]"] = round(v_tot, 3)
                h = heading if port == "run" else seg.get("branch_dir")
                out["nodes"][sid]["ports"][port] = disp_P(st_p["P"])
                for c in children.get((sid, port), []):
                    visit(c, st_p, seg.get("size", size), h, p_path)
            return
        else:
            raise ValueError(f"Unknown segment type {t}")

        out["nodes"][sid]["P_out"] = disp_P(st2["P"])
        out["nodes"][sid]["T_out"] = disp_T(st2["T"])
        for c in children.get((sid, "out"), []):
            visit(c, st2, seg.get("size", size), new_h, path)

    st0 = {"P": P0, "T": T0, "pr": pr0}
    path0 = [{"x": 0.0, "P": P0, "T": T0, "E": 0.0, "N": 0.0,
              "z": float(inlet.get("elev", 0) or 0), "seg": "INLET"}]
    visit(root, st0, inlet["size"], inlet.get("heading", "east"), path0)

    for r in rows:
        for k, v in r.items():
            if isinstance(v, float) and ("dP" in k or k.startswith("P ") or k.startswith("T ")):
                r[k] = round(v, 4)

    out["summary"] = {
        "fluid": fluid.name, "phase": fluid.phase, "P_unit": Pu, "dP_unit": dPu, "T_unit": Tu,
        "len_unit": lu, "P_in": float(inlet["P"]), "T_in": float(inlet["T"]), "m_kgh": m_in * 3600,
        "rho_in": pr0["rho"], "Q_m3h_in": m_in / pr0["rho"] * 3600,
        "inlet_label": inlet.get("label") or "Inlet", "inlet_elev": float(inlet.get("elev", 0) or 0),
        "coolprop": HAS_COOLPROP, "feasible": not flags["infeasible"],
    }
    out["warnings"] = list(dict.fromkeys(out["warnings"]))
    out["ok"] = True
    return out


# ==========================================================================
# 8. CATALOGUE for the UI
# ==========================================================================
def catalog():
    fluids = {n: f.phase for n, f in builtin_fluids().items()}
    fluids["Custom liquid"] = "liquid"
    fluids["Custom gas"] = "gas"
    if HAS_COOLPROP:
        for n in COOLPROP_FLUIDS:
            fluids["CoolProp: " + n] = "auto"
    sizes = {}
    for nps, od in PIPE_OD.items():
        sizes[nps] = {sch: round(od - 2 * WALL[sch][nps], 2) for sch in WALL}
    return {
        "fluids": fluids, "pressure_units": list(PRESSURE_UNITS), "temp_units": TEMP_UNITS,
        "flow_units": list(FLOW_UNITS), "length_units": list(LENGTH_UNITS),
        "nps": list(PIPE_OD), "schedules": list(WALL), "id_table": sizes,
        "materials": MATERIALS,
        "elbows": {k: v[0] for k, v in ELBOWS.items()},
        "fittings": {k: v[0] for k, v in FITTINGS.items()},
        "thermal": thermal.catalog(), "coolprop": HAS_COOLPROP, "directions": list(DIRS),
    }
