"""Environment delta-beta analysis: mass matching, bootstrap, permutation.

Implements the statistics fixed in report/preregistration_tng300.md
(Addenda 1-4):

* delta_beta_global  : beta of ALL satellites inside R200 pooled over the
  hosts of a half (primary statistic).
* delta_beta_profile : mean over 8 equal-count radial bins of the stacked
  beta(r) profile (legacy TNG100 definition).
* Mass matching within slices of log10 M200 (default 0.1 dex).
* Host-level bootstrap (matching repeated in every resample) and a
  within-mass-slice permutation null for the primary statistic.
"""

import numpy as np

from kinematics import radial_tangential_split


# ----------------------------------------------------------------------
# Per-host preparation
# ----------------------------------------------------------------------
def prepare_hosts(halos):
    """Precompute per-host quantities used by all later statistics.

    Parameters
    ----------
    halos : list of dict
        Each dict has "pos", "vel", "r200", "logm" (as in satellite_data).

    Returns
    -------
    dict
        logm, n, and per-host sufficient statistics for the global beta
        (s1 = sum(vr - shift), s2 = sum((vr - shift)^2), t2 = sum(vt^2)),
        plus per-host arrays r (scaled by R200), vr, vt for the profile.
    """
    r_list, vr_list, vt_list = [], [], []
    for h in halos:
        r, vr, vt = radial_tangential_split(h["pos"], h["vel"])
        r_list.append(r / h["r200"])
        vr_list.append(vr)
        vt_list.append(vt)

    shift = float(np.mean(np.concatenate(vr_list)))  # variance is shift-invariant
    return {
        "logm": np.array([h["logm"] for h in halos], dtype=float),
        "n": np.array([len(v) for v in vr_list], dtype=float),
        "s1": np.array([np.sum(v - shift) for v in vr_list]),
        "s2": np.array([np.sum((v - shift) ** 2) for v in vr_list]),
        "t2": np.array([np.sum(v ** 2) for v in vt_list]),
        "r": r_list,
        "vr": vr_list,
        "vt": vt_list,
    }


# ----------------------------------------------------------------------
# Statistics
# ----------------------------------------------------------------------
def global_beta(hosts, idx):
    """Pooled beta of all satellites of the hosts in ``idx``.

    Uses the second moment for tangential dispersion and the variance with
    ddof=1 around the pooled mean for radial dispersion. ``idx`` may
    contain repeated indices (bootstrap).

    Parameters
    ----------
    hosts : dict
        Output of prepare_hosts.
    idx : ndarray of int

    Returns
    -------
    float
        beta, or NaN if fewer than 3 satellites or zero radial variance.
    """
    n = hosts["n"][idx].sum()
    if n < 3:
        return np.nan
    s1 = hosts["s1"][idx].sum()
    s2 = hosts["s2"][idx].sum()
    sr2 = (s2 - s1 ** 2 / n) / (n - 1.0)
    st2 = hosts["t2"][idx].sum() / n
    if sr2 <= 0:
        return np.nan
    return 1.0 - st2 / (2.0 * sr2)


def profile_betas(r, vr, vt, n_bins=8, n_min=20):
    """Equal-count radial beta profile (fast, sort-based).

    Reproduces the legacy ``beta_from_rvrvt`` of the TNG100 notebook:
    quantile edges, bins [edge_i, edge_i+1), NaN where n < n_min.

    Parameters
    ----------
    r, vr, vt : ndarray
        Scaled radius, radial velocity, tangential speed of the pooled
        satellites.
    n_bins : int
    n_min : int

    Returns
    -------
    ndarray, shape (n_bins,)
        beta per radial bin (NaN where undefined).
    """
    order = np.argsort(r, kind="stable")
    r_s = r[order]
    vr_s = vr[order] - np.mean(vr)
    vt2_s = vt[order] ** 2

    edges = np.quantile(r_s, np.linspace(0.0, 1.0, n_bins + 1))
    edges[-1] += 1e-10
    b = np.searchsorted(r_s, edges, side="left")
    b[0] = 0
    b[-1] = len(r_s)

    c1 = np.concatenate([[0.0], np.cumsum(vr_s)])
    c2 = np.concatenate([[0.0], np.cumsum(vr_s ** 2)])
    ct = np.concatenate([[0.0], np.cumsum(vt2_s)])

    out = np.full(n_bins, np.nan)
    for i in range(n_bins):
        lo, hi = b[i], b[i + 1]
        n = hi - lo
        if n < n_min:
            continue
        s1 = c1[hi] - c1[lo]
        s2 = c2[hi] - c2[lo]
        sr2 = (s2 - s1 ** 2 / n) / (n - 1.0)
        if sr2 <= 0:
            continue
        st2 = (ct[hi] - ct[lo]) / n
        out[i] = 1.0 - st2 / (2.0 * sr2)
    return out


def profile_mean_beta(hosts, idx, n_bins=8, n_min=20):
    """Mean over radial bins of the stacked beta profile of hosts ``idx``."""
    r = np.concatenate([hosts["r"][i] for i in idx])
    vr = np.concatenate([hosts["vr"][i] for i in idx])
    vt = np.concatenate([hosts["vt"][i] for i in idx])
    b = profile_betas(r, vr, vt, n_bins=n_bins, n_min=n_min)
    return np.nan if np.all(np.isnan(b)) else float(np.nanmean(b))


# ----------------------------------------------------------------------
# Mass matching
# ----------------------------------------------------------------------
def slice_ids(logm, slice_width=0.1):
    """Integer mass-slice label for each host (slices anchored at 0)."""
    return np.floor(np.asarray(logm) / slice_width).astype(int)


def match_by_mass(logm, high, rng, slice_width=0.1):
    """Mass-match the high- and low-environment halves.

    Within each slice of log10 M200, keep equal numbers of high and low
    hosts by randomly dropping the surplus of the larger half.

    Parameters
    ----------
    logm : ndarray
    high : ndarray of bool
        True for the high-environment half.
    rng : numpy.random.Generator
    slice_width : float
        Slice width in dex.

    Returns
    -------
    hi_idx, lo_idx : ndarray of int
        Indices of the kept hosts in each half (equal length).
    """
    logm = np.asarray(logm)
    high = np.asarray(high, dtype=bool)
    sid = slice_ids(logm, slice_width)
    hi_keep, lo_keep = [], []
    for s in np.unique(sid):
        in_s = np.where(sid == s)[0]
        h = in_s[high[in_s]]
        l = in_s[~high[in_s]]
        k = min(len(h), len(l))
        if k == 0:
            continue
        hi_keep.append(rng.choice(h, size=k, replace=False))
        lo_keep.append(rng.choice(l, size=k, replace=False))
    if not hi_keep:
        return np.array([], dtype=int), np.array([], dtype=int)
    return np.concatenate(hi_keep), np.concatenate(lo_keep)


# ----------------------------------------------------------------------
# Delta-beta for one set of hosts
# ----------------------------------------------------------------------
def delta_beta_pair(hosts, hi_idx, lo_idx, with_profile=True):
    """delta_beta (high - low) for the given host index sets.

    Returns
    -------
    dict
        global, profile (NaN if with_profile is False).
    """
    out = {"global": global_beta(hosts, hi_idx) - global_beta(hosts, lo_idx),
           "profile": np.nan}
    if with_profile:
        out["profile"] = (profile_mean_beta(hosts, hi_idx)
                          - profile_mean_beta(hosts, lo_idx))
    return out


# ----------------------------------------------------------------------
# Bootstrap
# ----------------------------------------------------------------------
def bootstrap_delta_beta(hosts, high, matched, n_boot=1000, slice_width=0.1,
                         seed=42):
    """Host-level bootstrap of delta_beta.

    Hosts are resampled with replacement from the WHOLE bin; each host
    keeps its environment label (the median split is not recomputed), and
    mass matching, if requested, is repeated inside every resample.

    Parameters
    ----------
    hosts : dict
        Output of prepare_hosts.
    high : ndarray of bool
    matched : bool
    n_boot : int
    slice_width : float
    seed : int

    Returns
    -------
    ndarray, shape (n_boot, 2)
        Columns: delta_beta_global, delta_beta_profile.
    """
    rng = np.random.default_rng(seed)
    n = len(high)
    out = np.full((n_boot, 2), np.nan)
    for b in range(n_boot):
        draw = rng.integers(0, n, size=n)
        if matched:
            hi_l, lo_l = match_by_mass(hosts["logm"][draw], high[draw], rng,
                                       slice_width)
            hi_idx, lo_idx = draw[hi_l], draw[lo_l]
        else:
            hi_idx, lo_idx = draw[high[draw]], draw[~high[draw]]
        if len(hi_idx) == 0 or len(lo_idx) == 0:
            continue
        res = delta_beta_pair(hosts, hi_idx, lo_idx)
        out[b] = (res["global"], res["profile"])
    return out


# ----------------------------------------------------------------------
# Permutation null (primary statistic)
# ----------------------------------------------------------------------
def permutation_null(hosts, high, n_perm=10000, slice_width=0.1, seed=42):
    """Null distribution of the matched delta_beta_global.

    Environment labels are shuffled within each mass slice (which keeps
    the mass distribution and the number of high/low hosts per slice),
    then matched and delta_beta_global recomputed.

    Returns
    -------
    ndarray, shape (n_perm,)
    """
    rng = np.random.default_rng(seed)
    high = np.asarray(high, dtype=bool)
    sid = slice_ids(hosts["logm"], slice_width)
    groups = []
    for s in np.unique(sid):
        in_s = np.where(sid == s)[0]
        nh = int(high[in_s].sum())
        nl = len(in_s) - nh
        if min(nh, nl) > 0:
            groups.append((in_s, nh, min(nh, nl)))

    out = np.full(n_perm, np.nan)
    for p in range(n_perm):
        hi_parts, lo_parts = [], []
        for in_s, nh, k in groups:
            perm = rng.permutation(in_s)
            hi_parts.append(perm[:k])          # random k of the 'high' label
            lo_parts.append(perm[nh:nh + k])   # random k of the 'low' label
        hi_idx = np.concatenate(hi_parts)
        lo_idx = np.concatenate(lo_parts)
        out[p] = global_beta(hosts, hi_idx) - global_beta(hosts, lo_idx)
    return out


def permutation_pvalues(observed, null):
    """One-sided permutation p-values with the +1 correction.

    Returns
    -------
    p_greater, p_less : float
        P(null >= observed) and P(null <= observed).
    """
    null = null[~np.isnan(null)]
    n = len(null)
    p_greater = (np.sum(null >= observed) + 1.0) / (n + 1.0)
    p_less = (np.sum(null <= observed) + 1.0) / (n + 1.0)
    return float(p_greater), float(p_less)


def holm_correction(pvals):
    """Holm step-down adjusted p-values (same order as input)."""
    p = np.asarray(pvals, dtype=float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * p[i])
        adj[i] = min(1.0, running)
    return adj


# ----------------------------------------------------------------------
# Full analysis of one mass bin
# ----------------------------------------------------------------------
def analyse_mass_bin(halos, high, slice_width=0.1, n_boot=1000, n_perm=10000,
                     min_per_half=50, seed=42):
    """Run the pre-registered delta_beta analysis for one mass bin.

    Parameters
    ----------
    halos : list of dict
        Hosts of this mass bin (pos, vel, r200, logm).
    high : ndarray of bool
        High-environment flag per host (Addendum 2 median split).
    slice_width : float
    n_boot, n_perm : int
    min_per_half : int
        Testability threshold, applied after matching.
    seed : int

    Returns
    -------
    dict
        n_hosts_low/high (unmatched), n_matched_per_half, testable, and for
        'matched' and 'unmatched': global/profile estimate, bootstrap sigma
        and 68% interval; plus permutation p-values for the matched global
        statistic (if testable).
    """
    high = np.asarray(high, dtype=bool)
    hosts = prepare_hosts(halos)
    rng = np.random.default_rng(seed)

    hi_m, lo_m = match_by_mass(hosts["logm"], high, rng, slice_width)
    res = {
        "n_hosts_low": int((~high).sum()),
        "n_hosts_high": int(high.sum()),
        "n_matched_per_half": int(len(hi_m)),
        "mean_logm_low_matched": float(hosts["logm"][lo_m].mean()) if len(lo_m) else np.nan,
        "mean_logm_high_matched": float(hosts["logm"][hi_m].mean()) if len(hi_m) else np.nan,
    }
    res["testable"] = bool(len(hi_m) >= min_per_half)

    variants = {
        "unmatched": (np.where(high)[0], np.where(~high)[0]),
        "matched": (hi_m, lo_m),
    }
    for name, (hi_idx, lo_idx) in variants.items():
        if len(hi_idx) == 0 or len(lo_idx) == 0:
            res[name] = None
            continue
        est = delta_beta_pair(hosts, hi_idx, lo_idx)
        boot = bootstrap_delta_beta(hosts, high, matched=(name == "matched"),
                                    n_boot=n_boot, slice_width=slice_width,
                                    seed=seed + (1 if name == "matched" else 2))
        entry = {}
        for j, stat in enumerate(["global", "profile"]):
            col = boot[:, j]
            col = col[~np.isnan(col)]
            entry[stat] = {
                "delta_beta": float(est[stat]),
                "sigma_boot": float(np.std(col, ddof=1)) if len(col) > 1 else np.nan,
                "ci68": (float(np.percentile(col, 16)), float(np.percentile(col, 84)))
                if len(col) > 1 else (np.nan, np.nan),
                "n_boot_valid": int(len(col)),
            }
        res[name] = entry

    if res["testable"] and res["matched"] is not None:
        null = permutation_null(hosts, high, n_perm=n_perm,
                                slice_width=slice_width, seed=seed + 3)
        p_gt, p_lt = permutation_pvalues(res["matched"]["global"]["delta_beta"], null)
        res["perm"] = {"p_greater": p_gt, "p_less": p_lt,
                       "null_sigma": float(np.nanstd(null)), "n_perm": int(n_perm)}
    else:
        res["perm"] = None
    return res


# ----------------------------------------------------------------------
# Toy data for validation
# ----------------------------------------------------------------------
def make_toy_hosts(n_hosts=600, delta_true=0.0, mass_slope=1.0, base_beta=0.0,
                   env_mass_corr=0.6, seed=0):
    """Synthetic hosts with a KNOWN environment effect on beta.

    beta_i = base_beta + mass_slope * (logm_i - 3.4) + delta_true * high_i.
    Environment correlates with mass (so an unmatched split is biased) but
    the matched split should recover delta_true.

    Returns
    -------
    halos : list of dict (pos, vel, r200, logm)
    high : ndarray of bool
    """
    rng = np.random.default_rng(seed)
    logm = rng.uniform(3.0, 3.9, n_hosts)
    z = (logm - logm.mean()) / logm.std()
    env_latent = env_mass_corr * z + np.sqrt(1 - env_mass_corr ** 2) * rng.normal(size=n_hosts)
    high = env_latent > np.median(env_latent)

    halos = []
    for i in range(n_hosts):
        n_sat = 8 + rng.poisson(12)
        r = rng.uniform(0.05 ** 3, 1.0, n_sat) ** (1.0 / 3.0)
        cos_t = rng.uniform(-1, 1, n_sat)
        sin_t = np.sqrt(1 - cos_t ** 2)
        phi = rng.uniform(0, 2 * np.pi, n_sat)
        pos = np.column_stack([r * sin_t * np.cos(phi),
                               r * sin_t * np.sin(phi), r * cos_t])
        r_hat = pos / np.linalg.norm(pos, axis=1, keepdims=True)
        a = rng.normal(size=(n_sat, 3))
        a -= (a * r_hat).sum(axis=1, keepdims=True) * r_hat
        a_hat = a / np.linalg.norm(a, axis=1, keepdims=True)
        b_hat = np.cross(r_hat, a_hat)

        beta_i = base_beta + mass_slope * (logm[i] - 3.4) + delta_true * float(high[i])
        beta_i = min(beta_i, 0.9)
        sig_t = np.sqrt(2.0 * (1.0 - beta_i))
        vel = (rng.normal(0, 1.0, n_sat)[:, None] * r_hat
               + rng.normal(0, sig_t / np.sqrt(2), n_sat)[:, None] * a_hat
               + rng.normal(0, sig_t / np.sqrt(2), n_sat)[:, None] * b_hat)
        halos.append({"pos": pos, "vel": vel, "r200": 1.0, "logm": float(logm[i])})
    return halos, high
