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


---
# Addendum 1 (2026-10-03): satellite mass floor and sample rules

Written after the TNG300-1 group catalog was downloaded, and before any
TNG300 satellite kinematics or beta values were computed.

- Satellite floor: SubhaloMass >= 0.2 (code units of 1e10 Msun/h, i.e.
  2e9 Msun/h). This is ~50 DM-particle masses in TNG300-1 (m_DM = 0.00398)
  and ~400 in TNG100-1 (m_DM = 5.056e-4).
- The SAME floor is applied to TNG100 and TNG300. TNG100 is re-analysed
  with this floor so the two simulations are compared like-for-like.
- Host halos are selected as before (mass bins above, GroupNsubs >= 6),
  but a host is kept only if it has >= 5 satellites inside R200 AFTER the
  floor is applied.
- All other definitions and the decision rule in the original
  pre-registration are unchanged.


---
# Addendum 2 (2026-10-03): testability rule and median-split definition

Written after the TNG300 satellite catalog was built (floor 0.2, >= 5
satellites; counts per bin known) and BEFORE any beta or delta_beta was
computed in TNG300.

- Observed host counts after the floor: bin 0 = 17, bin 1 = 237,
  bin 2 = 3404, bin 3 = 10382, bin 4 = 2507.
- A mass bin is TESTABLE only if each environment half contains >= 50
  hosts (comparable to the ~62-65 hosts per half in the TNG100 cluster bin
  that motivated this test). Untestable bins are reported as untestable,
  not as passes or failures. Bin 0 is expected to be untestable.
- Holm correction for the secondary prediction (delta_beta < 0) is applied
  over the testable bins among bins 0-3 only.
- Median split: hosts with environment count above the median are 'high',
  below are 'low'; hosts exactly at the median are assigned by a seeded
  random draw (seed = 42) so the halves are as equal as possible. The same
  rule is used for the TNG100 re-analysis.
- Hosts above the top mass edge (15 hosts with log10 M200 >= 4.5) are
  excluded from the primary test and reported separately.
- The low-mass bins are a selected sample (hosts with >= 5 resolved
  satellites above the floor) and are not representative of all halos in
  those mass ranges. This limitation will be stated in the report.


---
# Addendum 3 (2026-10-03): mass matching and definition of delta_beta

Written BEFORE any beta or delta_beta was computed in TNG300.

Motivation (diagnostic of host mass only, no beta): in the TNG300 cluster
bin the high-environment half is more massive than the low half (mean
log10 M200 = 3.470 vs 3.287; satellite-weighted 3.816 vs 3.440; Spearman
between environment count and log M200 = +0.31). Bins 1-3 show gaps of
0.02-0.04 dex. Because beta increases with host mass, an unmatched delta_beta
can partly reflect mass.

Definitions:
- delta_beta = beta_global(high) - beta_global(low), where beta_global is
  computed from all satellites inside R200, pooled over hosts in the half
  (second moment for tangential, variance with ddof=1 for radial).
  Radial beta profiles are shown for illustration, not tested.
- PRIMARY statistic: mass-matched delta_beta. Within each 0.1-dex slice of
  log10 M200, equal numbers of high- and low-environment hosts are kept by
  randomly dropping surplus hosts of the larger half (seed 42). Applied to
  all testable bins. The unmatched delta_beta is reported as a sensitivity
  check; if the two disagree, the matched result decides.
- Testability (>= 50 hosts per half) is evaluated after matching.
- Uncertainty: 1000 bootstrap resamples of hosts (with replacement, whole
  bin), repeating the matching in each resample.
- Permutation null: 10,000 shuffles of environment labels within each mass
  slice, re-matching and recomputing delta_beta each time.
- Robustness (only if primary passes): slice width 0.05 dex.


---
# Addendum 4 (2026-10-03): legacy TNG100 definition of delta_beta

Written BEFORE any beta or delta_beta was computed in TNG300.

Finding: inspection of the TNG100 notebook (environment_split_within_mass_bins)
shows the exploratory TNG100 delta_beta was the difference of the unweighted
means over the 8 equal-count radial bins of the stacked beta(r) profile
(np.nanmean of per-radial-bin beta, n_min = 20), with the median split
assigning hosts exactly at the median to the LOW half. No uncertainty on
delta_beta was computed.

Changes / additions:
- The PRIMARY test is unchanged: mass-matched delta_beta_global (Addendum 3).
- For every test (each testable bin; matched and unmatched) a second
  statistic, delta_beta_profile = mean over the 8 radial bins of stacked
  beta(r) in the high half minus the same in the low half, is also reported
  with bootstrap errors. This is the TNG100 legacy definition and allows
  direct comparison with the original numbers.
- If delta_beta_global and delta_beta_profile have opposite signs in bin 4,
  the result is reported as definition-dependent.
- TNG100 re-analysis: first reproduce the original TNG100 numbers from the
  original catalog (no floor, ties to low half) as a code check; then
  re-run with the Addendum 1-3 rules (floor, matching, balanced split) and
  both statistics.
