# TNG Anisotropy Atlas

Velocity anisotropy of satellite galaxy populations in IllustrisTNG:
a reproducible measurement pipeline and population atlas across halo
mass, environment, and galaxy-population bins.

**Status:** Complete.  
See `notebooks/01_toy_validation.ipynb` for
pipeline validation on synthetic halos with known anisotropy.
and `02_real_data.ipynb` for illustrisTNG data based analysis. 

## Goal

Measure satellite velocity anisotropy β(r) = 1 − σ_t²/2σ_r² in
IllustrisTNG-100 at z = 0, validated against literature benchmarks
(Wojtak et al. 2013; standard N-body halo profiles), as a baseline for
testing alternative dynamical frameworks.

## Structure

- `notebooks/` — reproducible analysis workflows
- `src/` — modular pipeline code (kinematics, statistics, plotting)
- `figures/` — final figures
- `report/` — write-up (LaTeX)

## Data

IllustrisTNG public data release (Nelson et al. 2019), TNG100-1,
snapshot 99, accessed via `illustris_python`.

## License

MIT
