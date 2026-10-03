# Pre-registration: TNG300 replication of the environment sign reversal

Date: 2026-10-03 (written before any TNG300 satellite data were analysed)

## Hypothesis (from TNG100, exploratory)
At fixed host mass, the effect of environment on satellite velocity
anisotropy beta changes sign with host mass: delta_beta < 0 (denser
environment -> more tangential) in bins 0-3, delta_beta > 0 in the
cluster bin (bin 4). TNG100 cluster bin had only ~62-65 halos per half.

## Definitions (fixed now)
- delta_beta = beta(high-env) - beta(low-env), median split within each
  mass bin, from the stacked R200-scaled profile (equal-count radial bins).
- Mass bins: same physical edges as TNG100,
  log10(M200 / 1e10 Msun/h) = (1.0,1.3), (1.3,1.7), (1.7,2.2), (2.2,3.0), (3.0,4.5).
- Environment: neighbour count within 5 Mpc/h, neighbour threshold
  log10(M200) >= 1.0.
- Errors: halo-level bootstrap (resampling host halos), 1000 resamples.
- Satellite subhalo mass floor: a single physical-mass floor resolved in
  both TNG100 and TNG300, chosen from resolution considerations BEFORE
  looking at any delta_beta result, and applied to both simulations.

## Primary test (single confirmatory test)
Delta_beta in bin 4 is > 0 at > 2 sigma in TNG300.
Secondary: delta_beta < 0 in bins 0-3 (Holm-corrected).

## Decision rule
- Primary + secondary pass, and robustness checks pass: robust within the
  TNG family; consider a preprint.
- Primary passes but a robustness check fails: report as definition-dependent.
- Primary fails: report as not reproduced in TNG300.

## Robustness checks (run only if the primary passes)
Environment radius 2/5/10 Mpc/h; neighbour thresholds; finer mass sub-bins
and mass-matched comparison in bin 4; tercile and continuous-Spearman
variants; satellite-count and inner-radius cuts; TNG300 octant split for
cosmic variance; optional TNG300-Dark.

## Framing
Standard-gravity baseline, 3D velocities. Any positive result is
"suggestive, worth further investigation", not a discovery.
