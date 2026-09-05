"""TNG data I/O: authenticated download and catalog loading."""

import os
from pathlib import Path

import numpy as np
import requests


def get_api_key():
    """Retrieve TNG API key from Colab Secrets or environment variables.

    Returns
    -------
    str or None
    """
    try:
        from google.colab import userdata
        key = userdata.get("TNG_API_KEY")
        if key:
            return key
    except Exception:
        pass
    return os.environ.get("TNG_API_KEY")


def _download_file(url, out_path, api_key, chunk_size=8192):
    """Download a single file with API-Key header."""
    headers = {"API-Key": api_key}
    r = requests.get(url, headers=headers, stream=True, timeout=60)
    r.raise_for_status()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wb") as f:
        for chunk in r.iter_content(chunk_size=chunk_size):
            if chunk:
                f.write(chunk)
    return out_path


def download_groupcat(snapshot, base_path, api_key=None, max_chunks=128):
    """Download all fof_subhalo_tab chunks for a snapshot.

    Parameters
    ----------
    snapshot : int
        Snapshot number (e.g. 99).
    base_path : str or Path
        Local root where ``groups_099/`` will be created.
    api_key : str, optional
        TNG API key. If None, reads from Colab Secrets or env.
    max_chunks : int
        Safety limit on chunk enumeration.

    Returns
    -------
    list of Path
        Paths to downloaded files.
    """
    if api_key is None:
        api_key = get_api_key()
    if api_key is None:
        raise ValueError(
            "TNG API key not found. Set TNG_API_KEY in Colab Secrets "
            "or as an environment variable."
        )

    base_url = (
        f"http://www.tng-project.org/api/TNG100-1/files/groupcat-{snapshot}/"
    )
    out_dir = Path(base_path) / f"groups_{snapshot:03d}"
    out_dir.mkdir(parents=True, exist_ok=True)

    downloaded = []
    for i in range(max_chunks):
        fname = f"fof_subhalo_tab_{snapshot:03d}.{i}.hdf5"
        url = base_url + fname
        out_file = out_dir / fname

        if out_file.exists():
            downloaded.append(out_file)
            continue

        try:
            _download_file(url, out_file, api_key)
            downloaded.append(out_file)
            print(f"[tng_io] downloaded {fname}")
        except requests.exceptions.HTTPError as exc:
            if exc.response.status_code == 404:
                print(f"[tng_io] stop at chunk {i} (404)")
                break
            raise

    return downloaded


def load_catalogs(snapshot, base_path, group_fields=None, sub_fields=None):
    """Load group and subhalo catalogs via illustris_python.

    Parameters
    ----------
    snapshot : int
    base_path : str or Path
    group_fields : list of str, optional
    sub_fields : list of str, optional

    Returns
    -------
    group_cat, sub_cat, header : dict, dict, dict
    """
    import illustris_python as il

    base_path = str(base_path)
    group_cat = il.groupcat.loadHalos(base_path, snapshot, fields=group_fields)
    sub_cat = il.groupcat.loadSubhalos(base_path, snapshot, fields=sub_fields)
    header = il.groupcat.loadHeader(base_path, snapshot)
    return group_cat, sub_cat, header
