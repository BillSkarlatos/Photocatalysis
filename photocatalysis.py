#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quantum-Enhanced Photocatalysis — Equal Dose, Correct Calibration
Outputs:
  - Four separate figures (PNG, 300 dpi):
      1) rate_vs_flux.png
      2) yield_vs_flux.png
      3) enhancement_vs_flux.png
      4) bleaching_vs_flux.png   (includes: bleaching fraction, bqy, damage per product)
  - CSV with per-flux statistics: results.csv
"""

# Headless backend (prevents QSocketNotifier issues)
import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from enum import Enum
import pandas as pd


# ----------------------------
#   States & Parameters
# ----------------------------
class State(Enum):
    G = 0   # ground
    S1 = 1  # singlet
    T1 = 2  # triplet
    X = 3   # bleached / deactivated


@dataclass
class Params:
    # Optics
    beam_area_m2: float = 1e-13   # ~180 nm radius (tight, realistic for confocal)

    # Kinetics (Ru(bpy)3^2+-like exemplar; adjust if you have measurements)
    tau_s1_s: float = 5e-9
    tau_t1_s: float = 6e-7
    phi_isc: float = 0.6
    phi_react: float = 0.4

    # Classical 2PA cross-section (per-molecule rate r2 = sigma2 * I^2)
    sigma2_class_m4s: float = 1e-53   # ~1000 GM

    # Entangled calibration (see k_ent below)
    E_ref: float = 30.0               # target enhancement at reference flux
    Phi_ref: float = 3e12             # reference photon flux (photons/s)
    Phi_sat: float = 2e13             # entanglement roll-off scale (photons/s)

    # Bleaching per 2PA excitation
    p_bleach_2pa: float = 5e-6


# ----------------------------
#   Simulator (tau-leaping)
# ----------------------------
class Sim:
    def __init__(self, n_catalysts=200_000, params: Params = Params(), seed=7):
        self.n = int(n_catalysts)
        self.p = params
        self.rng = np.random.default_rng(seed)
        self.k_ent = self._calibrate_k_ent()

    def set_rng(self, rng: np.random.Generator):
        self.rng = rng

    def _calibrate_k_ent(self) -> float:
        """
        Solve for k_ent so that at Phi_ref:
            r2_ent(Phi_ref) = E_ref * r2_class(Phi_ref)
        Classical: r2_class = sigma2 * (Phi/A)^2
        Entangled: r2_ent  = k_ent * (Phi/2)
        => k_ent = 2 * E_ref * sigma2 * Phi_ref / A^2
        """
        p = self.p
        return 2.0 * p.E_ref * p.sigma2_class_m4s * p.Phi_ref / (p.beam_area_m2 ** 2)

    def r2_class(self, flux: float) -> float:
        I = flux / self.p.beam_area_m2
        return self.p.sigma2_class_m4s * (I ** 2)

    def r2_ent(self, flux: float) -> float:
        pair_flux = flux / 2.0
        roll = 1.0 / (1.0 + flux / self.p.Phi_sat)  # roll-off at high flux
        return self.k_ent * pair_flux * roll

    def _reset(self):
        self.state = np.full(self.n, State.G.value, dtype=np.uint8)
        self.products = 0
        self.absorbed = 0
        self.bleached = 0

    def _pois(self, lam):
        return self.rng.poisson(lam) if lam > 0 else 0

    def _choose(self, mask, k):
        idx = np.flatnonzero(mask)
        if k <= 0 or idx.size == 0:
            return np.empty(0, dtype=int)
        if k >= idx.size:
            return idx
        return self.rng.choice(idx, size=k, replace=False)

    def simulate(self, flux: float, duration: float, mode: str, n_steps: int = 600) -> dict:
        """Tau-leaping with integrated probabilities; mode in {'classical','entangled'}."""
        self._reset()
        dt = duration / n_steps

        r2 = self.r2_ent(flux) if mode == "entangled" else self.r2_class(flux)
        p_exc2 = 1.0 - np.exp(-r2 * dt)

        p_s1 = 1.0 - np.exp(-dt / self.p.tau_s1_s)
        p_t1 = 1.0 - np.exp(-dt / self.p.tau_t1_s)

        for _ in range(n_steps):
            if not np.any(self.state != State.X.value):
                break

            # Ground -> S1 via 2PA
            gmask = (self.state == State.G.value)
            Ng = int(gmask.sum())
            if Ng:
                n2 = self._pois(Ng * p_exc2)
                idx = self._choose(gmask, n2)
                if idx.size:
                    self.absorbed += 2 * idx.size
                    # bleaching on excitation
                    bl_mask = self.rng.random(idx.size) < self.p.p_bleach_2pa
                    if bl_mask.any():
                        bl = idx[bl_mask]
                        self.state[bl] = State.X.value
                        self.bleached += bl.size
                    survivors = idx[self.state[idx] == State.G.value]
                    if survivors.size:
                        self.state[survivors] = State.S1.value

            # S1 decay & branching
            s1mask = (self.state == State.S1.value)
            N1 = int(s1mask.sum())
            if N1:
                ndec = self._pois(N1 * p_s1)
                dec = self._choose(s1mask, ndec)
                if dec.size:
                    u = self.rng.random(dec.size)
                    f = (u < 0.10)  # fluorescence 10%
                    isc = (u >= 0.10) & (u < 0.10 + self.p.phi_isc)
                    rest = ~(f | isc)
                    if f.any():   self.state[dec[f]] = State.G.value
                    if isc.any(): self.state[dec[isc]] = State.T1.value
                    if rest.any():self.state[dec[rest]] = State.G.value

            # T1 decay: reaction vs non-productive
            t1mask = (self.state == State.T1.value)
            Nt = int(t1mask.sum())
            if Nt:
                ntd = self._pois(Nt * p_t1)
                td = self._choose(t1mask, ntd)
                if td.size:
                    u = self.rng.random(td.size)
                    react = td[u < self.p.phi_react]
                    if react.size:
                        self.products += react.size
                    self.state[td] = State.G.value

        eps = 1e-30
        return {
            "flux": flux,
            "duration": duration,
            "rate": self.products / max(duration, eps),
            "yield_abs": self.products / max(self.absorbed, 1),
            "eff_inc": self.products / max(flux * duration, eps),
            "bleach_frac": self.bleached / self.n,
            "bleaches": self.bleached,
            "absorptions": self.absorbed,
            "products": self.products,
            # expected classical 2PA events per run (for reliability mask)
            "exp_class_events": self.n * (1.0 - np.exp(-self.r2_class(flux) * duration)),
        }


# ----------------------------
#   Paired runs, bootstrap
# ----------------------------
def paired_compare(sim: Sim, flux: float, duration: float, seed: int):
    base = np.random.default_rng(seed)
    sC = int(base.integers(0, 2**63 - 1)); sE = int(base.integers(0, 2**63 - 1))
    sim.set_rng(np.random.default_rng(sC)); rc = sim.simulate(flux, duration, "classical")
    sim.set_rng(np.random.default_rng(sE)); re = sim.simulate(flux, duration, "entangled")
    return rc, re

def mean_ci_boot(x, alpha=0.05, B=2000, rng=None):
    rng = np.random.default_rng() if rng is None else rng
    x = np.asarray(x, float)
    m = x.mean()
    if x.size < 2:
        return float(m), (float(m), float(m))
    boots = np.empty(B)
    n = x.size
    for b in range(B):
        idx = rng.integers(0, n, n)
        boots[b] = x[idx].mean()
    lo, hi = np.quantile(boots, [alpha/2, 1-alpha/2])
    return float(m), (float(lo), float(hi))

def ratio_ci_boot(num, den, alpha=0.05, B=4000, rng=None):
    rng = np.random.default_rng() if rng is None else rng
    num = np.asarray(num, float); den = np.asarray(den, float)
    eps = 1e-30
    logs = np.log((num + eps) / (den + eps))
    m = float(np.exp(logs.mean()))
    if logs.size < 2:
        return m, (m, m)
    boots = np.empty(B); n = logs.size
    for b in range(B):
        idx = rng.integers(0, n, n)
        boots[b] = np.exp(logs[idx].mean())
    lo, hi = np.quantile(boots, [alpha/2, 1-alpha/2])
    return float(m), (float(lo), float(hi))

def binom_ci_approx(k, n, alpha=0.05):
    if n <= 0: return (0.0, 0.0)
    p = k / n
    se = np.sqrt(max(p*(1-p), 1e-30) / n)
    return (max(0.0, p - 1.96*se), min(1.0, p + 1.96*se))


# ----------------------------
#   Sweep (equal incident dose)
# ----------------------------
def sweep(sim: Sim, fluxes, duration=1.0, trials=32, seed0=20251005, reliable_min=50.0):
    out_rows = []  # for CSV (one row per flux with summary stats)

    summary = {
        "flux": np.array(fluxes, float),
        # means
        "c_rate": [], "e_rate": [], "c_yld": [], "e_yld": [],
        "c_eff": [], "e_eff": [],
        "enh": [],
        # CIs
        "c_rate_ci": [], "e_rate_ci": [],
        "c_yld_ci": [], "e_yld_ci": [],
        "c_eff_ci": [], "e_eff_ci": [],
        "enh_ci": [],
        # bleaching (aggregated over trials to kill rare-event noise)
        "c_bleach_frac": [], "e_bleach_frac": [],
        "c_bleach_frac_ci": [], "e_bleach_frac_ci": [],
        "bqy_c": [], "bqy_e": [],            # bleaches per absorbed photon
        "bqy_c_ci": [], "bqy_e_ci": [],
        "dpp_c": [], "dpp_e": [],            # damage per product (bleaches/product)
        "reliable": [],
    }

    for i, F in enumerate(fluxes):
        # collect paired trial results
        rate_c = []; rate_e = []
        yld_c = [];  yld_e = []
        eff_c = [];  eff_e = []
        exp_class = []

        # aggregate counts across trials (for stable bleaching stats)
        Bc = 0; Be = 0   # bleaches
        Ac = 0; Ae = 0   # absorptions
        Pc = 0; Pe = 0   # products

        for t in range(trials):
            rc, re = paired_compare(sim, F, duration, seed0 + 1000*i + t)
            rate_c.append(rc["rate"]);  rate_e.append(re["rate"])
            yld_c.append(rc["yield_abs"]); yld_e.append(re["yield_abs"])
            eff_c.append(rc["eff_inc"]);  eff_e.append(re["eff_inc"])
            exp_class.append(rc["exp_class_events"])

            Bc += rc["bleaches"]; Be += re["bleaches"]
            Ac += rc["absorptions"]; Ae += re["absorptions"]
            Pc += rc["products"];    Pe += re["products"]

        # means + CIs
        m, ci = mean_ci_boot(rate_c); summary["c_rate"].append(m); summary["c_rate_ci"].append(ci)
        m, ci = mean_ci_boot(rate_e); summary["e_rate"].append(m); summary["e_rate_ci"].append(ci)
        m, ci = mean_ci_boot(yld_c);  summary["c_yld"].append(m);  summary["c_yld_ci"].append(ci)
        m, ci = mean_ci_boot(yld_e);  summary["e_yld"].append(m);  summary["e_yld_ci"].append(ci)
        m, ci = mean_ci_boot(eff_c);  summary["c_eff"].append(m);  summary["c_eff_ci"].append(ci)
        m, ci = mean_ci_boot(eff_e);  summary["e_eff"].append(m);  summary["e_eff_ci"].append(ci)

        # enhancement on paired ratios
        enh, enh_ci = ratio_ci_boot(rate_e, rate_c)
        summary["enh"].append(enh); summary["enh_ci"].append(enh_ci)

        # reliability mask (expected classical excitations per run)
        reliable = (np.mean(exp_class) >= reliable_min)
        summary["reliable"].append(reliable)

        # bleaching fraction (aggregate across trials)
        frac_c = Bc / max(sim.n * trials, 1)
        frac_e = Be / max(sim.n * trials, 1)
        # CI using binomial approx on N*trials trials
        ci_c = binom_ci_approx(Bc, sim.n * trials)
        ci_e = binom_ci_approx(Be, sim.n * trials)
        summary["c_bleach_frac"].append(frac_c); summary["c_bleach_frac_ci"].append(ci_c)
        summary["e_bleach_frac"].append(frac_e); summary["e_bleach_frac_ci"].append(ci_e)

        # Bleaching quantum yield (per absorbed photon) + CI
        bqy_c = Bc / max(Ac, 1); bqy_e = Be / max(Ae, 1)
        summary["bqy_c"].append(bqy_c); summary["bqy_e"].append(bqy_e)
        summary["bqy_c_ci"].append(binom_ci_approx(Bc, Ac)); summary["bqy_e_ci"].append(binom_ci_approx(Be, Ae))

        # Damage per product (bleaches per product)
        dpp_c = Bc / max(Pc, 1e-30); dpp_e = Be / max(Pe, 1e-30)
        summary["dpp_c"].append(dpp_c); summary["dpp_e"].append(dpp_e)

        # ---- CSV row (flatten CIs) ----
        out_rows.append({
            "flux": F,
            "rate_class_mean": summary["c_rate"][-1],   "rate_class_lo": ci_or(summary["c_rate_ci"][-1], 0), "rate_class_hi": ci_or(summary["c_rate_ci"][-1], 1),
            "rate_ent_mean":   summary["e_rate"][-1],   "rate_ent_lo":   ci_or(summary["e_rate_ci"][-1], 0), "rate_ent_hi":   ci_or(summary["e_rate_ci"][-1], 1),
            "yield_class_mean":summary["c_yld"][-1],    "yield_class_lo":ci_or(summary["c_yld_ci"][-1], 0),  "yield_class_hi":ci_or(summary["c_yld_ci"][-1], 1),
            "yield_ent_mean":  summary["e_yld"][-1],    "yield_ent_lo":  ci_or(summary["e_yld_ci"][-1], 0),  "yield_ent_hi":  ci_or(summary["e_yld_ci"][-1], 1),
            "eff_class_mean":  summary["c_eff"][-1],    "eff_class_lo":  ci_or(summary["c_eff_ci"][-1], 0),  "eff_class_hi":  ci_or(summary["c_eff_ci"][-1], 1),
            "eff_ent_mean":    summary["e_eff"][-1],    "eff_ent_lo":    ci_or(summary["e_eff_ci"][-1], 0),  "eff_ent_hi":    ci_or(summary["e_eff_ci"][-1], 1),
            "enh_mean":        summary["enh"][-1],      "enh_lo":        ci_or(summary["enh_ci"][-1], 0),     "enh_hi":        ci_or(summary["enh_ci"][-1], 1),
            "reliable": int(reliable),
            "bleach_frac_class": frac_c, "bleach_frac_class_lo": ci_c[0], "bleach_frac_class_hi": ci_c[1],
            "bleach_frac_ent":   frac_e, "bleach_frac_ent_lo":   ci_e[0], "bleach_frac_ent_hi":   ci_e[1],
            "bqy_class": bqy_c, "bqy_class_lo": summary["bqy_c_ci"][-1][0], "bqy_class_hi": summary["bqy_c_ci"][-1][1],
            "bqy_ent":   bqy_e, "bqy_ent_lo":   summary["bqy_e_ci"][-1][0], "bqy_ent_hi":   summary["bqy_e_ci"][-1][1],
            "damage_per_product_class": dpp_c,
            "damage_per_product_ent":   dpp_e,
            "sum_bleaches_class": Bc, "sum_absorptions_class": Ac, "sum_products_class": Pc,
            "sum_bleaches_ent":   Be, "sum_absorptions_ent":   Ae, "sum_products_ent":   Pe,
        })

    # convert lists to arrays for plotting convenience
    for k in list(summary.keys()):
        if k == "flux": continue
        summary[k] = np.array(summary[k], dtype=object if k.endswith("_ci") else float)

    # write CSV
    df = pd.DataFrame(out_rows)
    df.to_csv("results.csv", index=False)

    return summary


def ci_or(ci_tuple, i):
    """Helper: extract CI bound when stored as (lo,hi) tuple."""
    if isinstance(ci_tuple, (list, tuple)) and len(ci_tuple) == 2:
        return float(ci_tuple[i])
    return float(ci_tuple)


# ----------------------------
#   Plotting (four files)
# ----------------------------
def band(ax, x, ci, color=None):
    lo = np.array([t[0] for t in ci], float)
    hi = np.array([t[1] for t in ci], float)
    ax.fill_between(x, lo, hi, alpha=0.25, linewidth=0, color=color)

def plot_rate(summary, fname="rate_vs_flux.png"):
    F = summary["flux"]; eps = 1e-30
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    band(ax, F, summary["c_rate_ci"]); band(ax, F, summary["e_rate_ci"])
    ax.plot(F, summary["c_rate"]+eps, "o-", label="Classical")
    ax.plot(F, summary["e_rate"]+eps, "s-", label="Entangled")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Photon flux (photons s$^{-1}$)")
    ax.set_ylabel("Reaction rate (products s$^{-1}$)")
    ax.set_title("Reaction Rate vs Flux (equal incident photons)")
    ax.grid(alpha=0.3, which="both"); ax.legend()
    fig.savefig(fname, dpi=300)

def plot_yield(summary, fname="yield_vs_flux.png"):
    F = summary["flux"]; eps = 1e-30
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    band(ax, F, summary["c_yld_ci"]); band(ax, F, summary["e_yld_ci"])
    ax.plot(F, summary["c_yld"]+eps, "o-", label="Classical")
    ax.plot(F, summary["e_yld"]+eps, "s-", label="Entangled")
    ax.set_xscale("log")
    ax.set_xlabel("Photon flux (photons s$^{-1}$)")
    ax.set_ylabel("Quantum yield (products / absorbed photon)")
    ax.set_title("Quantum Yield vs Flux")
    ax.grid(alpha=0.3, which="both"); ax.legend()
    fig.savefig(fname, dpi=300)

def plot_enh(summary, fname="enhancement_vs_flux.png"):
    F = summary["flux"]; eps = 1e-30
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    mask = np.array(summary["reliable"], bool)
    # reliable
    if mask.any():
        enh_ci = summary["enh_ci"][mask]
        lo = np.array([c[0] for c in enh_ci]); hi = np.array([c[1] for c in enh_ci])
        ax.fill_between(F[mask], lo, hi, color="green", alpha=0.2, linewidth=0)
        ax.plot(F[mask]+eps, summary["enh"][mask]+eps, "^-", color="green", label="reliable")
    # unreliable (hollow markers, not connected)
    if (~mask).any():
        ax.scatter(F[~mask]+eps, summary["enh"][~mask]+eps, facecolors="none", edgecolors="green", label="unreliable")
    ax.axhline(1.0, ls="--", c="k", lw=1)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Photon flux (photons s$^{-1}$)")
    ax.set_ylabel("Enhancement (Entangled / Classical)")
    ax.set_title("Quantum Enhancement (95% CI; hollow = unreliable)")
    ax.grid(alpha=0.3, which="both"); ax.legend()
    fig.savefig(fname, dpi=300)

def plot_bleach(summary, fname="bleaching_vs_flux.png"):
    F = summary["flux"]; eps = 1e-30
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 10.0), constrained_layout=True)

    # A) Bleached fraction (aggregated)
    ax = axes[0]
    band(ax, F, summary["c_bleach_frac_ci"]); band(ax, F, summary["e_bleach_frac_ci"])
    ax.plot(F, summary["c_bleach_frac"]+eps, "o-", label="Classical")
    ax.plot(F, summary["e_bleach_frac"]+eps, "s-", label="Entangled")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylabel("Bleached fraction")
    ax.set_title("Photobleaching vs Flux (aggregated across trials)")
    ax.grid(alpha=0.3, which="both"); ax.legend()

    # B) Bleaching quantum yield (bleaches per absorbed photon)
    ax = axes[1]
    band(ax, F, summary["bqy_c_ci"]); band(ax, F, summary["bqy_e_ci"])
    ax.plot(F, summary["bqy_c"]+eps, "o-", label="Classical")
    ax.plot(F, summary["bqy_e"]+eps, "s-", label="Entangled")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_ylabel("Bleaches / absorbed photon")
    ax.set_title("Bleaching Quantum Yield")
    ax.grid(alpha=0.3, which="both"); ax.legend()

    # C) Damage per product
    ax = axes[2]
    ax.plot(F, summary["dpp_c"]+eps, "o-", label="Classical")
    ax.plot(F, summary["dpp_e"]+eps, "s-", label="Entangled")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Photon flux (photons s$^{-1}$)")
    ax.set_ylabel("Bleaches / product")
    ax.set_title("Photodamage per Product (lower is better)")
    ax.grid(alpha=0.3, which="both"); ax.legend()

    fig.savefig(fname, dpi=300)


# ----------------------------
#   Main
# ----------------------------
def main():
    # Flux range: covers reliable low-flux region while avoiding ultra-rare events
    flux = np.logspace(11, 13.8, 14)   # ~1e11 .. 6e13 photons/s
    duration = 1.0                      # equal incident photons per flux point
    trials = 32

    sim = Sim(n_catalysts=200_000, params=Params(), seed=7)
    summary = sweep(sim, flux, duration=duration, trials=trials, seed0=20251005, reliable_min=50.0)

    # Four independent figures
    plot_rate(summary, "rate_vs_flux.png")
    plot_yield(summary, "yield_vs_flux.png")
    plot_enh(summary, "enhancement_vs_flux.png")
    plot_bleach(summary, "bleaching_vs_flux.png")

    # Brief console summary for the reliable region
    mask = np.array(summary["reliable"], bool)
    if mask.any():
        i = int(np.argmax(summary["enh"][mask]))
        idx = np.flatnonzero(mask)[i]
        print("\n=== SUMMARY (reliable region) ===")
        print(f"Peak enhancement: {summary['enh'][idx]:.2f}× at flux {summary['flux'][idx]:.2e} photons/s")
    else:
        print("\n[Note] No reliable points met the threshold; increase duration/molecules.")
    print("Saved figures: rate_vs_flux.png, yield_vs_flux.png, enhancement_vs_flux.png, bleaching_vs_flux.png")
    print("Saved CSV: results.csv")

if __name__ == "__main__":
    main()
