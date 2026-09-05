"""Host-halo selection and satellite catalog construction."""

import numpy as np


def periodic_distance(pos_a, pos_b, box_size):
    """Compute periodic separation vector pos_a - pos_b in a cubic box.

    Parameters
    ----------
    pos_a, pos_b : ndarray, shape (..., 3)
    box_size : float
        Box side length in same units as positions.

    Returns
    -------
    delta : ndarray
        Separation vector with periodic wrapping applied.
    """
    delta = pos_a - pos_b
    delta = np.where(delta > box_size / 2.0, delta - box_size, delta)
    delta = np.where(delta < -box_size / 2.0, delta + box_size, delta)
    return delta


def select_host_halos(group_cat, logm_edges, min_nsubs=6):
    """Select host halos in mass bins with enough subhalos.

    Parameters
    ----------
    group_cat : dict
        Must contain Group_M_Crit200, Group_R_Crit200, GroupPos,
        GroupVel, GroupNsubs, GroupFirstSub.
    logm_edges : array-like, shape (n_bins, 2)
        log10(M / (1e10 Msun/h)) bin edges.
    min_nsubs : int
        Minimum GroupNsubs required (central + satellites).

    Returns
    -------
    list of dict
        One entry per mass bin with keys: indices, mass, r200,
        pos, vel, first_sub, logm.
    """
    mass = group_cat["Group_M_Crit200"]
    logm = np.log10(mass)
    nsubs = group_cat["GroupNsubs"]

    hosts = []
    for low, high in logm_edges:
        mask = (logm >= low) & (logm < high) & (nsubs >= min_nsubs)
        idx = np.where(mask)[0]

        hosts.append(
            {
                "indices": idx,
                "mass": mass[idx],
                "r200": group_cat["Group_R_Crit200"][idx],
                "pos": group_cat["GroupPos"][idx],
                "vel": group_cat["GroupVel"][idx],
                "first_sub": group_cat["GroupFirstSub"][idx],
                "logm": logm[idx],
            }
        )
    return hosts


def build_satellite_catalog(sub_cat, host_idx, host_props, box_size,
                            min_stellar_mass=0.0):
    """Build satellite catalog for a single host halo.

    Parameters
    ----------
    sub_cat : dict
        Subhalo catalog with required fields.
    host_idx : int
        FOF group index.
    host_props : dict
        Single-host properties: pos, vel, r200, first_sub.
    box_size : float
        Periodic box size.
    min_stellar_mass : float
        Minimum stellar mass in code units (10^10 Msun/h).
        Default 0.0 keeps all flagged subhalos.

    Returns
    -------
    rel_pos : ndarray, shape (N_sat, 3)
    rel_vel : ndarray, shape (N_sat, 3)
    props : dict
        Satellite properties: r, stellar_mass, sfr, flag.
    """
    grnr = sub_cat["SubhaloGrNr"]
    in_group = grnr == host_idx

    # Exclude central
    central = host_props["first_sub"]
    is_sat = in_group & (np.arange(len(grnr)) != central)

    # Quality & stellar-mass cuts
    is_sat &= sub_cat["SubhaloFlag"] > 0
    stellar = sub_cat["SubhaloMassType"][:, 4]
    is_sat &= stellar >= min_stellar_mass

    if not np.any(is_sat):
        return None, None, None

    sat_pos = sub_cat["SubhaloPos"][is_sat]
    sat_vel = sub_cat["SubhaloVel"][is_sat]

    # Periodic position offset
    rel_pos = periodic_distance(sat_pos, host_props["pos"], box_size)
    r = np.sqrt(np.sum(rel_pos ** 2, axis=1))

    # Velocity offset (no periodic wrapping for velocities)
    rel_vel = sat_vel - host_props["vel"]

    # R200 cut
    within = r < host_props["r200"]
    if not np.any(within):
        return None, None, None

    props = {
        "r": r[within],
        "stellar_mass": stellar[is_sat][within],
        "sfr": sub_cat["SubhaloSFR"][is_sat][within],
        "flag": sub_cat["SubhaloFlag"][is_sat][within],
    }

    return rel_pos[within], rel_vel[within], props
