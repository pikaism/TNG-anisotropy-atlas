"""TNG data I/O: authenticated download and catalog loading."""

import os
import time
from pathlib import Path

import numpy as np
import requests

HDF5_MAGIC = b"\x89HDF\r\n\x1a\n"


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
    except Exception as e:
        print(f"[tng_io] Colab secrets lookup failed: {type(e).__name__}: {e}")
    return os.environ.get("TNG_API_KEY")


def _is_valid_hdf5(path, min_size=1024):
    """Check the downloaded file has the correct HDF5 magic bytes."""
    try:
        if path.stat().st_size < min_size:
            return False
        with open(path, "rb") as f:
            return f.read(8) == HDF5_MAGIC
    except OSError:
        return False


def _download_file(url, out_path, api_key, chunk_size=8192,
                    max_retries=6, backoff=15, max_backoff=90):
    """Download a single file with api-key header, retries, and HDF5 validation."""
    headers = {"api-key": api_key}
    url = url.replace("http://", "https://")

    for attempt in range(max_retries):
        try:
            r = requests.get(url, headers=headers, stream=True, timeout=60)
            if r.status_code != 200:
                raise RuntimeError(f"status {r.status_code}")

            out_path.parent.mkdir(parents=True, exist_ok=True)
            with open(out_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=chunk_size):
                    if chunk:
                        f.write(chunk)

            if not _is_valid_hdf5(out_path):
                out_path.unlink(missing_ok=True)
                raise RuntimeError("downloaded content is not valid HDF5")

            return out_path

        except (RuntimeError, requests.exceptions.RequestException) as e:
            if attempt < max_retries - 1:
                wait = min(backoff * (attempt + 1), max_backoff)
                print(f"  [retry {attempt+1}/{max_retries}] {out_path.name}: "
                      f"{e} -- waiting {wait}s")
                time.sleep(wait)
            else:
                raise


def download_groupcat(snapshot, base_path, api_key=None):
    """Download all groupcat chunks for a snapshot using the TNG API file listing.

    Parameters
    ----------
    snapshot : int
        Snapshot number (e.g. 99).
    base_path : str or Path
        Local root where ``groups_099/`` will be created.
    api_key : str, optional
        TNG API key. If None, reads from Colab Secrets or env.

    Returns
    -------
    list of Path
        Paths to downloaded (or already-present, verified) files.
    """
    if api_key is None:
        api_key = get_api_key()
    if api_key is None:
        raise ValueError(
            "TNG API key not found. Set TNG_API_KEY in Colab Secrets "
            "or as an environment variable."
        )

    headers = {"api-key": api_key}
    listing_url = f"https://www.tng-project.org/api/TNG100-1/files/groupcat-{snapshot}/"
    r = requests.get(listing_url, headers=headers, timeout=30)
    r.raise_for_status()
    file_urls = r.json()["files"]

    out_dir = Path(base_path) / f"groups_{snapshot:03d}"
    out_dir.mkdir(parents=True, exist_ok=True)

    downloaded = []
    for url in file_urls:
        fname = url.split("/")[-1]
        out_file = out_dir / fname

        if out_file.exists() and _is_valid_hdf5(out_file):
            downloaded.append(out_file)
            continue

        _download_file(url, out_file, api_key)
        downloaded.append(out_file)
        print(f"[tng_io] downloaded {fname}")

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
