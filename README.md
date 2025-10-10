# Quantum‑Enhanced Photocatalysis — Equal Dose, Correct Calibration

**Script:** `photocatalysis.py`
**Purpose:** Simulate and compare classical two‑photon absorption (2PA) photocatalysis against an *entangled‑photon–enabled* pathway under **equal incident photon dose**. The simulator sweeps photon flux, runs paired Monte‑Carlo trials, and produces publication‑ready figures and a CSV of per‑flux statistics.

---

## What this code does

* Implements a minimal four‑state photocatalyst model: **G** (ground), **S1** (singlet), **T1** (triplet), **X** (bleached/deactivated).
* Models **classical 2PA**: rate per molecule ( r^\text{class} = \sigma_2 I^2 ).
* Models an **entangled pathway**: ( r^\text{ent} = k_\text{ent},(\Phi/2),R(\Phi) ) with a tunable roll‑off ( R(\Phi)=1/(1+\Phi/\Phi_\text{sat}) ).
* **Calibrates** (k_\text{ent}) so that at a reference flux (\Phi_\text{ref}), the entangled rate matches a target enhancement (E_\text{ref}) over classical.
* Uses **tau‑leaping** Monte‑Carlo with paired RNG seeds to reduce variance when taking **ratios** (enhancement).
* Tracks reaction products, absorbed photons, and **photobleaching** events; computes confidence intervals via bootstrap.
* Sweeps flux (equal exposure time → equal incident photons per point) and writes **four figures** + **CSV**.

---

## Outputs

On completion the script writes:

1. `rate_vs_flux.png` — Reaction rate vs photon flux (log–log), with 95% CIs.
2. `yield_vs_flux.png` — Quantum yield (products/absorbed photon) vs flux, with 95% CIs.
3. `enhancement_vs_flux.png` — Ratio **(Entangled / Classical)** vs flux. Filled line = reliable region; hollow markers = low‑count region.
4. `bleaching_vs_flux.png` — Three stacked panels:

   * Bleached fraction (aggregated across trials) with 95% CIs.
   * Bleaching quantum yield (bleaches per absorbed photon) with 95% CIs.
   * Photodamage per product (bleaches/product) — lower is better.
5. `results.csv` — Per‑flux summary including means, CIs, reliability mask, and aggregated counts.

The console also prints a short summary of the **peak enhancement** within the reliability mask.

---

## Quick start

### Requirements

* Python ≥ 3.9
* Packages: `numpy`, `matplotlib`, `pandas`

```bash
python3 -m pip install numpy matplotlib pandas
```

### Run

```bash
python3 photocatalysis.py
```

Artifacts are saved in the working directory.

---

## How to interpret the model

* **Equal dose:** Each flux point uses the same **duration** (default 1 s), so incident photons (\propto \Phi \times t) are constant across the sweep.
* **Classical vs entangled:**

  * Classical 2PA scales (\propto I^2 = (\Phi/A)^2) with beam area (A).
  * Entangled term scales (\propto \Phi) at low flux and rolls off above `Phi_sat` to reflect practical degradation of entanglement advantages.
* **Calibration:** `k_ent` is set so that at `Phi_ref` the entangled event rate equals `E_ref ×` the classical rate. This enforces apples‑to‑apples comparison under your chosen reference.
* **States & branching:**

  * S1 decays by fluorescence (10%), ISC to T1 with probability `phi_isc`, or returns to G.
  * T1 reacts with probability `phi_react`; otherwise returns to G.
  * Bleaching occurs with probability `p_bleach_2pa` **at excitation**.

---

## Key parameters (edit in code via `Params`)

* **Optics**

  * `beam_area_m2` — Effective illuminated area. Smaller area → higher intensity for a given flux.
* **Kinetics** (Ru(bpy)₃²⁺‑like exemplar; adjust to taste)

  * `tau_s1_s`, `tau_t1_s` — Lifetimes of S1 and T1.
  * `phi_isc`, `phi_react` — ISC yield and reaction yield.
* **Classical 2PA**

  * `sigma2_class_m4s` — 2PA cross‑section (m⁴·s). ~1000 GM ≈ 10⁻⁵³ m⁴·s.
* **Entangled calibration**

  * `E_ref` — Target enhancement at `Phi_ref`.
  * `Phi_ref` — Reference flux (photons/s).
  * `Phi_sat` — Roll‑off scale where entangled advantage diminishes.
* **Photobleaching**

  * `p_bleach_2pa` — Probability of bleaching per 2PA excitation.

**Simulation controls** (edit calls to `sweep`/`main`):

* `n_catalysts` — Number of molecules simulated.
* `flux` — Log‑spaced array of flux values (see `main`).
* `duration` — Exposure time per flux point (equal dose).
* `trials` — Paired Monte‑Carlo replicates for CIs and ratios.
* `reliable_min` — Minimum expected classical excitations per run to mark a point as "reliable".

---

## Reproducibility & statistics

* **Paired RNG:** For each flux the classical and entangled runs are seeded from a common base to reduce ratio variance.
* **CIs:**

  * Means: non‑parametric bootstrap of the mean.
  * Ratios: bootstrap on log‑ratios (geometric mean & CI).
  * Bleaching fraction / bqy: normal approximation to binomial; aggregated over all trials to stabilize rare events.

---

## Modifying the sweep

In `main()`:

```python
flux = np.logspace(11, 13.8, 14)  # adjust range and density
duration = 1.0                     # equal dose per point
trials = 32                        # more trials → tighter CIs
```

You can also modify `reliable_min` (expected classical events per run) to shift the reliability mask.

---

## Output schema (`results.csv`)

Each row corresponds to one flux value and includes:

* `rate_*`, `yield_*`, `eff_*` — Means and 95% CI bounds for classical and entangled.
* `enh_*` — Enhancement (entangled/classical) with 95% CI.
* `reliable` — 1 if the classical expected events ≥ `reliable_min`.
* Bleaching metrics: fractions with CIs, **BQY** (bleaches per absorbed photon with CIs), and **damage per product**.
* Aggregated counts: total bleaches/absorptions/products across trials for each mode.

---

## Tips & gotchas

* Ensure `beam_area_m2` realistically matches your optical setup; rate scaling is sensitive to it.
* If **no reliable points** appear, increase `duration`, `n_catalysts`, or raise flux range.
* Use consistent units: flux in photons/s, area in m², lifetimes in s.
* For very low fluxes, events are rare; CIs widen and points may be flagged unreliable (hollow in the enhancement plot).

---

## Extending the model

* Replace the fixed fluorescence (10%) with a parameter, or make it wavelength‑dependent.
* Add oxygen quenching, triplet–triplet annihilation, or diffusion‑limited encounter kinetics.
* Allow bleaching during T1 chemistry, not only at excitation.
* Introduce spatial heterogeneity (hot spots) by sampling local intensity.

---

## License

Choose and add a license if you plan to share (e.g., MIT, BSD‑3‑Clause, Apache‑2.0).

---

## Citation

Pending . . .

---

## Acknowledgments

* This code path draws on standard tau‑leaping Monte‑Carlo ideas and common photophysics parameterizations; tune parameters to your specific chromophore/catalyst.
