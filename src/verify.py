"""
Delft University of Technology (TU Delft) hereby disclaims all copyright interest in the
program 'LocScale2' written by the Author(s).

Copyright (C) 2026 Alok Bharadwaj and Arjen J. Jakobi

The `locscale2 verify` command. Given an EMDB id from the LocScale-FEM paper, it opens the
published ChimeraX session (its FEM and pVDDT maps), fetches the EMDB map, and points the
LocScale-FEM tool at it so the user can re-run and compare.

Only the panels this bundle can reproduce are included: figures 6 and 7, and supplementary
figures 7 and 8. The hybrid-LocScale panels are deliberately excluded.
"""
import os
import re
import tempfile
import urllib.request

from chimerax.core.commands import CmdDesc, StringArg, quote_if_necessary
from chimerax.core.errors import UserError
from chimerax.map import Volume


# Public SURFdrive share holding the published sessions (read-only, no login).
SHARE_TOKEN = "MW8Bfdb3HD3P9Rq"
_DAV_URL = "https://surfdrive.surf.nl/public.php/dav/files/{token}/{name}"

# EBI half-map location (deposited under /other/ when the authors provided them).
_EBI_HALF = ("https://ftp.ebi.ac.uk/pub/databases/emdb/structures/"
             "EMD-{token}/other/emd_{token}_half_map_{n}.map.gz")

# EMDB id (exactly as it appears in the file name, zero-padded where the deposition is)
# -> published session file(s), figure panel first. See the share listing.
VERIFY_CATALOG = {
    #"13234": ["figure_6a_6b_13234.cxs"],
    #"17929": ["figure_6c_17929.cxs"],
    #"19999": ["figure_6d_6e_6f_19999.cxs"],
    #"0665":  ["figure_6g_0665.cxs"],
    "33888": ["figure_7a_33888.cxs"],
    "35193": ["figure_7b_35193.cxs"],
    "11231": ["figure_7c_11231.cxs", "supplementary_8d_11231.cxs"],
    "10366": ["figure_7c_10366.cxs", "supplementary_8e_10366.cxs"],
    "10577": ["figure_7c_10577.cxs", "supplementary_8c_10577.cxs"],
    #"15949": ["supplementary_7d_15949.cxs", "supplementary_7f_15949.cxs"],
    "7770":  ["supplementary_8b_7770.cxs"],
    "33394": ["supplementary_8a_33394.cxs"],
}


def _resolve_token(emdb):
    """Map user input ('13234', 'EMD-13234', '665') to a catalog key, or None."""
    digits = re.sub(r"\D", "", emdb or "")
    if not digits:
        return None
    target = int(digits)
    for token in VERIFY_CATALOG:
        if int(token) == target:
            return token
    return None


def _cache_path(name):
    cache_dir = os.path.join(tempfile.gettempdir(), "locscale2_verify")
    os.makedirs(cache_dir, exist_ok=True)
    return os.path.join(cache_dir, name)


def _download_share_file(session, name):
    """Download a file from the public share into a cache; reuse it if already present."""
    dest = _cache_path(name)
    if os.path.exists(dest) and os.path.getsize(dest) > 0:
        session.logger.info("Using cached {}.".format(name))
        return dest

    url = _DAV_URL.format(token=SHARE_TOKEN, name=urllib.request.quote(name))
    session.logger.status("Downloading {} ...".format(name))
    last = [0]

    def report(blocks, block_size, total):
        if total > 0:
            pct = min(100, int(100 * blocks * block_size / total))
            if pct >= last[0] + 10:
                last[0] = pct
                session.logger.status("Downloading {} ... {}%".format(name, pct))

    tmp = dest + ".part"
    try:
        urllib.request.urlretrieve(url, tmp, reporthook=report)
    except Exception as exc:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise UserError("Could not download {} from the share: {}".format(name, exc))
    os.replace(tmp, dest)
    session.logger.info("Downloaded {}.".format(name))
    return dest


def _emdb_symmetry(token):
    """Best-effort point group from the EMDB REST API; None on any problem."""
    try:
        import json
        url = "https://www.ebi.ac.uk/emdb/api/entry/EMD-{}".format(token)
        with urllib.request.urlopen(url, timeout=30) as resp:
            return _find_value(json.load(resp), "point_group")
    except Exception:
        return None


def _find_value(obj, key):
    """First string value stored under `key` anywhere in a nested dict/list, or None."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k == key and isinstance(v, str):
                return v
            found = _find_value(v, key)
            if found:
                return found
    elif isinstance(obj, list):
        for item in obj:
            found = _find_value(item, key)
            if found:
                return found
    return None


def _url_exists(url):
    """True if a HEAD request returns 200."""
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=30) as resp:
            return getattr(resp, "status", 200) == 200
    except Exception:
        return False


def _open_url_volume(session, url, name):
    from chimerax.core.commands import run
    opened = run(session, "open {}".format(url))
    vols = [m for m in (opened or []) if isinstance(m, Volume)]
    if not vols:
        return None
    vols[0].name = name
    return vols[0]


def _open_input_maps(session, token):
    """Open the two EMDB half maps (averaged input) when both exist, else the primary map.

    Returns ('halves', [v1, v2]) or ('single', [v]).
    """
    from chimerax.core.commands import run
    h1 = _EBI_HALF.format(token=token, n=1)
    h2 = _EBI_HALF.format(token=token, n=2)
    if _url_exists(h1) and _url_exists(h2):
        session.logger.status("Fetching half maps for EMD-{} ...".format(token))
        v1 = _open_url_volume(session, h1, "EMD-{} half map 1".format(token))
        v2 = _open_url_volume(session, h2, "EMD-{} half map 2".format(token))
        if v1 is not None and v2 is not None:
            return "halves", [v1, v2]
    session.logger.status("No half maps; fetching the primary EMD-{} map ...".format(token))
    opened = run(session, "open emdb:{}".format(token))
    vols = [m for m in (opened or []) if isinstance(m, Volume)]
    return "single", vols


def _ensure_tool(session):
    """The running LocScale-FEM tool, starting it if necessary."""
    from .tool import LocScale2Tool
    existing = session.tools.find_by_class(LocScale2Tool)
    if existing:
        return existing[0]
    from chimerax.core.commands import run
    run(session, "ui tool show LocScale2")
    existing = session.tools.find_by_class(LocScale2Tool)
    return existing[0] if existing else None


def _log_catalog(session):
    lines = ["LocScale-FEM verifiable EMDB entries (figures 6, 7; supplementary 7, 8):"]
    for token in sorted(VERIFY_CATALOG, key=lambda t: int(t)):
        lines.append("  EMD-{}  ({})".format(token, ", ".join(VERIFY_CATALOG[token])))
    session.logger.info("\n".join(lines))


def verify(session, emdb, panel=None):
    """Load an EMDB map and its published LocScale-FEM session for verification."""
    from chimerax.core.commands import run

    if emdb and emdb.strip().lower() in ("list", "?"):
        _log_catalog(session)
        return

    token = _resolve_token(emdb)
    if token is None:
        _log_catalog(session)
        raise UserError("EMDB '{}' is not in the verifiable set (figures 6, 7 and "
                        "supplementary 7, 8). See the list above.".format(emdb))

    sessions = VERIFY_CATALOG[token]
    if panel is not None:
        sessions = [s for s in sessions if panel in s] or sessions
    session_file = sessions[0]

    session.logger.info("LocScale-FEM verify: EMD-{}".format(token))

    # 1. Published session (the paper's FEM and pVDDT maps). Restoring a session replaces the
    #    current models, so it goes first and the tool is (re)shown afterwards.
    session.logger.warning("Opening the published session replaces the current ChimeraX scene.")
    path = _download_share_file(session, session_file)
    run(session, "open {}".format(quote_if_necessary(path)))
    session.logger.info("Opened published session {}.".format(session_file))
    if len(sessions) > 1:
        session.logger.info(
            "Other sessions for this map: {}. Open one with "
            "'locscale2 verify {} panel <name>'.".format(", ".join(sessions[1:]), token))

    # The maps from the restored session are the published references for the CC report.
    published = list(session.models.list(type=Volume))

    # 2. EMDB input map(s): half maps (averaged) when available, else the primary map.
    kind, vols = _open_input_maps(session, token)
    if not vols:
        raise UserError("Could not fetch EMD-{} maps from EMDB.".format(token))

    # 3. Show the tool, point it at the input, remember the references, fill in symmetry.
    tool = _ensure_tool(session)
    if tool is not None:
        tool._verify_reference = {"token": token, "published": published}
        if kind == "halves":
            tool._half1_menu.value = vols[0]
            tool._half2_menu.value = vols[1]
            session.logger.info("Loaded two half maps for EMD-{} (averaged at run time).".format(token))
        else:
            tool._map_menu.value = vols[0]
            session.logger.info("Loaded the primary EMD-{} map as the input.".format(token))
        pg = _emdb_symmetry(token)
        if pg and pg.upper() != "C1":
            tool._point_group_symmetry_menu.setText(pg)
            session.logger.info("Set point-group symmetry to {} (from EMDB).".format(pg))
    else:
        session.logger.warning("Could not open the LocScale-FEM tool; open it from the "
                               "Tools menu and select EMD-{} as the input.".format(token))

    session.logger.info(
        "Ready: EMD-{} is loaded in LocScale-FEM. Click 'Run feature enhancement' to "
        "reproduce; a cross-correlation with the published maps in {} is reported "
        "afterwards.".format(token, session_file))


verify_desc = CmdDesc(
    required=[("emdb", StringArg)],
    keyword=[("panel", StringArg)],
    synopsis="Load an EMDB map and its published LocScale-FEM session for verification")
