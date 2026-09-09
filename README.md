# TNG Anisotropy Atlas

A reproducible pipeline for measuring satellite galaxy velocity
anisotropy in IllustrisTNG, and an exploratory, validated analysis of
how it depends on halo mass, environment, and galaxy population.

**Status:** Complete.
See `notebooks/01_toy_validation.ipynb` for pipeline validation on
synthetic halos with known anisotropy, and `notebooks/02_real_data.ipynb`
for the IllustrisTNG-based analysis.

## Goal

Measure satellite velocity anisotropy β(r) = 1 − σ_t²/(2σ_r²) in
IllustrisTNG-100 at z = 0, validated against literature benchmarks
(Wojtak, Hansen & Hjorth 2013; standard N-body halo profiles). This is
a validated measurement pipeline and pilot analysis, not a claim of
novel physics — the value of this project is a reproducible pipeline,
honest validation, and a baseline for future comparisons with
alternative dynamical frameworks.

## Key findings

- Toy-model ensemble tests recover imposed β values to within ~1σ
  across β ∈ {-0.5, 0, 0.5}, confirming the pipeline is unbiased
  before use on real data.
- In the TNG100-1 sample, stacked β(r) shows a monotonic trend from
  tangentially-biased orbits at group scale (β ≈ -0.32 to -0.03) to
  radially-biased orbits at cluster scale (β ≈ +0.07 to +0.24),
  consistent with WHH13 and standard N-body expectations.
- Environment has a statistically significant effect on β independent
  of host-halo mass (partial Spearman ρ = -0.099, p = 4.19×10⁻⁷,
  N=2580 halos), with a sign reversal between group and cluster scale
  that we treat as suggestive rather than confirmed, given small
  sample size in the cluster bin.
- Quenched satellites are consistently more radially biased than
  star-forming ones, and quenched fraction declines monotonically
  from cluster center outward (94% → 70% in the cluster bin) —
  consistent with known ram-pressure/tidal quenching physics, used
  here as an independent validation of the pipeline.
- A halo mass function check confirms the analysis sample
  (log₁₀(M₂₀₀) ≥ 1.0) sits ~2.6 dex above the simulation's
  ~50-particle resolution limit.

## What this project contains

- **Toy-model validation**: synthetic halos with known input β,
  recovering the correct value via ensemble tests.
- **Real IllustrisTNG-100 data** (snapshot 99): host halos selected
  in 5 mass bins, satellite kinematics measured, stacked β(r) profiles
  with bootstrap errors (resampling host halos, not individual
  satellites, to respect intra-halo correlation).
- **Environment analysis**: neighbor-count-based environment measure,
  Spearman and partial-Spearman correlations with β, and a
  mass-controlled environment split.
- **Resolution checks**: a halo mass function confirming the analysis
  sample sits well above the resolution limit, and a documented
  stellar-mass resolution floor for resolved galaxies.
- **Galaxy-population splits**: full-subhalo vs. resolved-galaxy-only
  comparison, and star-forming/quenched splits, including a radial
  quenched-fraction profile in cluster-mass halos.

## Definitions

- **β(r)**: velocity anisotropy, 1 − σ_t²/(2σ_r²), where σ_t² is the
  second moment of tangential speed and σ_r² is the radial-velocity
  variance around the bin mean.
- **Radial normalization**: satellite distances are scaled by each
  host halo's own R200 (r/R200), so profiles from halos of different
  size can be stacked together.
- **Satellite**: any subhalo gravitationally bound to a host (via
  SubhaloGrNr), excluding the central subhalo, within R200 of the
  host center.
- **Host-halo mass bins**: five bins in log₁₀(M200 / 10¹⁰ M☉/h):
  [1.0, 1.3), [1.3, 1.7), [1.7, 2.2), [2.2, 3.0), [3.0, 4.5).
- **Star-forming / quenched**: classified via specific star formation
  rate (SFR / stellar mass), threshold 10⁻¹¹ yr⁻¹ (Franx et al. 2008).
- **Environment**: number of neighboring halos above
  log₁₀(M200) ≥ 1.0 within a 5 Mpc/h periodic search radius.
- **Partial Spearman correlation**: β vs. environment, controlling for
  host-halo mass via rank-residualization.

## Structure
src/
kinematics.py — radial_tangential_split(), beta_profile()
halo_sample.py — periodic_distance(), select_host_halos(), build_satellite_catalog()
tng_io.py — TNG data download and catalog loading
stats_utils.py — stacked/per-halo statistics, correlations, environment,
SF/quenched classification, radial quenched-fraction profile
plotting.py — reusable plotting functions
notebooks/
01_toy_validation.ipynb — pipeline validation on synthetic data
02_real_data.ipynb — full IllustrisTNG-100 analysis
figures/ — final figures (see below)
report/ — LaTeX write-up


## Figures

**Beta Environment Split**

![Beta Environment Split](figures/beta_environment_split.png)

**Beta Galaxy Vs Full**

![Beta Galaxy Vs Full](figures/beta_galaxy_vs_full.png)

**Beta Sf Quenched**

![Beta Sf Quenched](figures/beta_sf_quenched.png)

**Beta Vs Mass**

![Beta Vs Mass](figures/beta_vs_mass.png)

**Halo Mass Function**

![Halo Mass Function](figures/halo_mass_function.png)

**Quenched Fraction Radial**

![Quenched Fraction Radial](figures/quenched_fraction_radial.png)


## Reproducibility

1. Install dependencies:
```bash
   pip install numpy scipy matplotlib illustris_python
```
2. Set `TNG_API_KEY` as an environment variable or Colab Secret
   (obtain from the TNG public data release site).
3. Run `notebooks/01_toy_validation.ipynb` first (no data download
   required — synthetic data only).
4. Run `notebooks/02_real_data.ipynb`. This downloads ~450 groupcat
   files from the TNG API on first run (large; cached locally after).

Raw TNG HDF5 data files are downloaded, not committed (excluded via
`.gitignore` due to size). Built satellite catalogs are cached to
Google Drive in the notebook workflow but are not part of this repo.

Developed and run on Python 3.11 (Google Colab), numpy, scipy,
matplotlib, illustris_python.

## Data

IllustrisTNG public data release (Nelson et al. 2019), TNG100-1,
snapshot 99, accessed via `illustris_python` and the TNG public API.

## License

MIT
