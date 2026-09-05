"""Kinematic utilities for satellite anisotropy."""

import numpy as np


def radial_tangential_split(pos, vel):
    """Decompose velocities into radial and tangential components.

    Parameters
    ----------
    pos : ndarray, shape (N, 3)
        Satellite positions relative to host center.
    vel : ndarray, shape (N, 3)
        Satellite velocities relative to host center.

    Returns
    -------
    r : ndarray, shape (N,)
    vr : ndarray, shape (N,)
        Radial velocity (positive = outflow).
    vt : ndarray, shape (N,)
        Tangential speed magnitude.
    """
    r = np.sqrt(np.sum(pos ** 2, axis=1))
    r_hat = pos / r[:, None]
    vr = np.sum(vel * r_hat, axis=1)
    vt_vec = vel - vr[:, None] * r_hat
    vt = np.sqrt(np.sum(vt_vec ** 2, axis=1))
    return r, vr, vt


def beta_profile(pos, vel, n_bins=8, n_min=20):
    """Measure beta(r) in equal-count radial bins.

    beta = 1 - sigma_t^2 / (2 * sigma_r^2)

    Uses:
    - second moment <vt^2> for tangential dispersion
    - variance with ddof=1 around bin mean for radial dispersion

    Parameters
    ----------
    pos, vel : ndarray, shape (N, 3)
    n_bins : int
    n_min : int
        Minimum satellites per bin to report a value.

    Returns
    -------
    centers : ndarray
        Bin centers in radius.
    beta : ndarray
        beta per bin (NaN where n < n_min).
    beta_err : ndarray
        Delta-method uncertainties.
    edges : ndarray
        Radial bin edges.
    """
    r, vr, vt = radial_tangential_split(pos, vel)

    edges = np.quantile(r, np.linspace(0.0, 1.0, n_bins + 1))
    edges[-1] += 1e-10  # guard against float edge case

    centers = []
    beta_vals = []
    beta_errs = []

    for i in range(n_bins):
        mask = (r >= edges[i]) & (r < edges[i + 1])
        n = int(mask.sum())

        if n < n_min:
            centers.append(0.5 * (edges[i] + edges[i + 1]))
            beta_vals.append(np.nan)
            beta_errs.append(np.nan)
            continue

        vr_b = vr[mask]
        vt_b = vt[mask]

        # Tangential: second moment
        st2 = np.mean(vt_b ** 2)
        st2_err = st2 / np.sqrt(n)

        # Radial: variance around bin mean
        sr2 = np.var(vr_b, ddof=1)
        if sr2 == 0:
            centers.append(0.5 * (edges[i] + edges[i + 1]))
            beta_vals.append(np.nan)
            beta_errs.append(np.nan)
            continue
        sr2_err = sr2 * np.sqrt(2.0 / (n - 1))

        beta = 1.0 - st2 / (2.0 * sr2)

        # Delta-method error propagation
        var_beta = (
            (1.0 / (2.0 * sr2)) ** 2 * st2_err ** 2
            + (st2 / (2.0 * sr2 ** 2)) ** 2 * sr2_err ** 2
        )
        beta_err = np.sqrt(var_beta)

        centers.append(0.5 * (edges[i] + edges[i + 1]))
        beta_vals.append(beta)
        beta_errs.append(beta_err)

    return (
        np.array(centers),
        np.array(beta_vals),
        np.array(beta_errs),
        edges,
    )
