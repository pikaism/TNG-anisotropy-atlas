"""Stacked, per-halo, and correlation statistics for satellite anisotropy.

NOTE: filename intentionally keeps the "statisitics" typo (see repo
README/handoff notes) — do not silently rename it.
"""

import numpy as np
from scipy.stats import rankdata, pearsonr

from kinematics import radial_tangential_split


def stack_beta_profile(bin_halos, n_bins=8, n_min=20, n_bootstrap=200, seed=42):
    """Stack satellites across halos (scaled by R200) and compute beta(r)
    with bootstrap errors resampling over HOST HALOS (not individual
    satellites).

    Satellites within one halo are not independent draws, so
    bootstrapping individual satellites would understate the true
    uncertainty; resampling whole halos preserves that correlation
    structure.

    Parameters
    ----------
    bin_halos : list of dict
        Each dict has "pos", "vel", "r200" for one halo.
    n_bins : int
        Number of equal-count radial bins.
    n_min : int
        Minimum satellites required in a bin to report a value.
    n_bootstrap : int
        Number of halo-resampling bootstrap iterations.
    seed : int
        RNG seed for reproducibility.

    Returns
    -------
    centers : ndarray
    beta : ndarray
    beta_err : ndarray
        Bootstrap standard deviation per bin.
    edges : ndarray
    n_halos : int
    """
    rng = np.random.default_rng(seed)
    n_halos = len(bin_halos)

    def stacked_r_vr_vt(halo_subset):
        all_r, all_vr, all_vt = [], [], []
        for h in halo_subset:
            r, vr, vt = radial_tangential_split(h["pos"], h["vel"])
            all_r.append(r / h["r200"])
            all_vr.append(vr)
            all_vt.append(vt)
        return np.concatenate(all_r), np.concatenate(all_vr), np.concatenate(all_vt)

    def beta_from_rvrvt(r, vr, vt, n_bins, n_min):
        edges = np.quantile(r, np.linspace(0, 1, n_bins + 1))
        edges[-1] += 1e-10
        centers, betas = [], []
        for i in range(n_bins):
            mask = (r >= edges[i]) & (r < edges[i + 1])
            n = mask.sum()
            if n < n_min:
                centers.append(np.nan)
                betas.append(np.nan)
                continue
            vr_b, vt_b = vr[mask], vt[mask]
            st2 = np.mean(vt_b ** 2)
            sr2 = np.var(vr_b, ddof=1)
            if sr2 == 0:
                centers.append(np.nan)
                betas.append(np.nan)
                continue
            betas.append(1.0 - st2 / (2.0 * sr2))
            centers.append(0.5 * (edges[i] + edges[i + 1]))
        return np.array(centers), np.array(betas), edges

    r_all, vr_all, vt_all = stacked_r_vr_vt(bin_halos)
    centers, beta_main, edges = beta_from_rvrvt(r_all, vr_all, vt_all, n_bins, n_min)

    boot_betas = np.full((n_bootstrap, n_bins), np.nan)
    for b in range(n_bootstrap):
        idx = rng.integers(0, n_halos, size=n_halos)
        r_b, vr_b, vt_b = stacked_r_vr_vt([bin_halos[i] for i in idx])
        _, beta_b, _ = beta_from_rvrvt(r_b, vr_b, vt_b, n_bins, n_min)
        boot_betas[b] = beta_b

    return centers, beta_main, np.nanstd(boot_betas, axis=0), edges, n_halos


def compute_halo_level_beta(bin_halos, n_min=5):
    """Compute a single global (non-radially-binned) beta per halo.

    Used for per-halo correlation analysis against mass/environment, as
    opposed to stack_beta_profile's radially-resolved stacked profile.

    Parameters
    ----------
    bin_halos : list of dict
        Each dict has "pos", "vel", "r200", "logm", "host_idx".
    n_min : int
        Minimum satellites required to compute a per-halo beta.

    Returns
    -------
    list of dict
        One entry per halo with a valid beta: host_idx, logm, beta, n_sat.
    """
    results = []
    for h in bin_halos:
        r, vr, vt = radial_tangential_split(h["pos"], h["vel"])
        n = len(r)
        if n < n_min:
            continue
        sr2 = np.var(vr, ddof=1)
        st2 = np.mean(vt ** 2)
        if sr2 == 0:
            continue
        results.append({
            "host_idx": h["host_idx"],
            "logm": h["logm"],
            "beta": 1.0 - st2 / (2.0 * sr2),
            "n_sat": n,
        })
    return results


def partial_spearman(x, y, z):
    """Partial Spearman correlation between x and y, controlling for z,
    via rank-residualization.

    Ranks x, y, z; regresses ranked x and ranked y linearly on ranked z;
    Pearson-correlates the two residual series. Standard nonparametric
    approximation to a partial Spearman correlation.

    Parameters
    ----------
    x, y, z : array-like
        Equal-length 1D arrays.

    Returns
    -------
    rho : float
    p : float
        p-value from the Pearson test on the residuals (approximate;
        does not adjust degrees of freedom for the two regressions).
    """
    x, y, z = np.asarray(x), np.asarray(y), np.asarray(z)
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)

    def _residuals(a, b):
        A = np.column_stack([np.ones_like(b), b])
        coef, *_ = np.linalg.lstsq(A, a, rcond=None)
        return a - A @ coef

    rx_resid = _residuals(rx, rz)
    ry_resid = _residuals(ry, rz)
    return pearsonr(rx_resid, ry_resid)


def compute_environment_counts(group_cat, box_size, radius_mpc_h=5.0,
                                mass_threshold_logm=1.0):
    """Count neighboring halos within a fixed radius, using periodic wrapping.

    Parameters
    ----------
    group_cat : dict
    box_size : float, ckpc/h
    radius_mpc_h : float
        Search radius in Mpc/h (converted internally to ckpc/h).
    mass_threshold_logm : float
        Only count neighbors above this log10(M200) threshold.

    Returns
    -------
    ndarray, shape (N_groups,)
        Neighbor count per group (0 if isolated or below threshold).
    """
    from scipy.spatial import cKDTree

    pos = group_cat["GroupPos"]
    mass = group_cat["Group_M_Crit200"]
    logm = np.full_like(mass, -np.inf)
    valid_mass = mass > 0
    logm[valid_mass] = np.log10(mass[valid_mass])

    massive_mask = logm >= mass_threshold_logm
    massive_pos = pos[massive_mask]
    massive_idx = np.where(massive_mask)[0]

    radius_ckpc = radius_mpc_h * 1000.0
    tree = cKDTree(massive_pos, boxsize=box_size)

    counts = np.zeros(len(pos), dtype=int)
    neighbor_lists = tree.query_ball_point(massive_pos, r=radius_ckpc)
    for local_i, neighbors in enumerate(neighbor_lists):
        counts[massive_idx[local_i]] = len(neighbors) - 1
    return counts


def classify_sf_quenched(props, ssfr_threshold=1e-11):
    """Classify satellites as star-forming or quenched via specific SFR.

    Parameters
    ----------
    props : dict
        Must contain "sfr" and "stellar_mass" (code units: 1e10 Msun/h).
    ssfr_threshold : float
        sSFR threshold in 1/yr; below this is "quenched". 1e-11/yr is a
        commonly used literature threshold (e.g. Franx et al. 2008).

    Returns
    -------
    is_quenched : ndarray of bool
    """
    stellar_mass_msun = props["stellar_mass"] * 1e10
    with np.errstate(divide="ignore", invalid="ignore"):
        ssfr = props["sfr"] / stellar_mass_msun
    return (ssfr < ssfr_threshold) | (props["sfr"] == 0)


def filter_to_galaxies(satellite_data, stellar_mass_floor):
    """Filter satellite_data to a resolved-galaxy stellar-mass floor.

    Drops halos left with fewer than 3 resolved-galaxy satellites.

    Parameters
    ----------
    satellite_data : list of dict
        Output of build_all_satellites (from halo_sample.py usage).
    stellar_mass_floor : float
        Minimum stellar mass in Msun (not code units).

    Returns
    -------
    list of dict
        Same structure, filtered.
    """
    filtered = []
    for bin_result in satellite_data:
        new_halos = []
        for h in bin_result["halos"]:
            stellar_mass_msun = h["props"]["stellar_mass"] * 1e10
            is_galaxy = stellar_mass_msun >= stellar_mass_floor
            if is_galaxy.sum() < 3:
                continue
            new_halos.append({
                **h,
                "pos": h["pos"][is_galaxy],
                "vel": h["vel"][is_galaxy],
                "props": {k: v[is_galaxy] for k, v in h["props"].items()},
            })
        filtered.append({"halos": new_halos})
    return filtered
