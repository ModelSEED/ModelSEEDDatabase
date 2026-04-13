import logging
import threading
import time

import requests

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
MAX_RETRIES = 4
NAME_BATCH_SIZE = 100  # names per POST request

logger = logging.getLogger("pubchem_validate")


def _extract_smiles(props):
    """Extract best SMILES from PubChem API response properties.

    PubChem returns SMILES under different keys depending on the endpoint:
    the key is often just 'SMILES' even when 'IsomericSMILES' was requested.
    Prefer IsomericSMILES (preserves stereochemistry) over CanonicalSMILES.
    """
    return (props.get("IsomericSMILES")
            or props.get("CanonicalSMILES")
            or props.get("SMILES", ""))


class RateLimiter:
    """Thread-safe token bucket rate limiter."""

    def __init__(self, rate=5.0):
        """rate: max requests per second."""
        self._interval = 1.0 / rate
        self._lock = threading.Lock()
        self._next_time = time.monotonic()

    def acquire(self):
        with self._lock:
            now = time.monotonic()
            if now < self._next_time:
                wait = self._next_time - now
                self._next_time += self._interval
            else:
                wait = 0
                self._next_time = now + self._interval
        if wait > 0:
            time.sleep(wait)


_rate_limiter = RateLimiter(rate=5.0)


def pubchem_request(method, url, **kwargs):
    """Make a PubChem request with retry and rate limiting."""
    for attempt in range(MAX_RETRIES):
        _rate_limiter.acquire()
        try:
            if method == "GET":
                resp = requests.get(url, timeout=30, **kwargs)
            else:
                resp = requests.post(url, timeout=30, **kwargs)
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code == 404:
                return None
            if resp.status_code in (429, 500, 502, 503, 504):
                wait = 2 ** (attempt + 1)
                logger.warning("PubChem %d, retrying in %ds...",
                               resp.status_code, wait)
                time.sleep(wait)
                continue
            return None
        except (requests.RequestException, ValueError):
            wait = 2 ** (attempt + 1)
            time.sleep(wait)
    return None


def query_xref(xref_id):
    """Query PubChem by xref ID. Returns (cid, inchikey, smiles) or None."""
    url = (f"{PUBCHEM_BASE}/compound/xref/RegistryID/"
           f"{xref_id}/property/InChIKey,IsomericSMILES/JSON")
    data = pubchem_request("GET", url)
    if data and "PropertyTable" in data:
        prop_list = data["PropertyTable"]["Properties"]
        if not prop_list:
            return None
        props = prop_list[0]
        return (
            str(props.get("CID", "")),
            props.get("InChIKey", ""),
            _extract_smiles(props),
        )
    return None


def query_inchikey(inchikey):
    """Query PubChem by InChIKey. Returns (cid, inchikey, smiles) or None."""
    url = (f"{PUBCHEM_BASE}/compound/inchikey/{inchikey}"
           f"/property/InChIKey,IsomericSMILES/JSON")
    data = pubchem_request("GET", url)
    if data and "PropertyTable" in data:
        prop_list = data["PropertyTable"]["Properties"]
        if not prop_list:
            return None
        props = prop_list[0]
        return (
            str(props.get("CID", "")),
            props.get("InChIKey", ""),
            _extract_smiles(props),
        )
    return None


def query_inchi(inchi_str):
    """Query PubChem by InChI string. Returns (cid, inchikey, smiles) or None."""
    url = (f"{PUBCHEM_BASE}/compound/inchi"
           f"/property/InChIKey,IsomericSMILES/JSON")
    data = pubchem_request("POST", url, data={"inchi": inchi_str})
    if data and "PropertyTable" in data:
        prop_list = data["PropertyTable"]["Properties"]
        if not prop_list:
            return None
        props = prop_list[0]
        return (
            str(props.get("CID", "")),
            props.get("InChIKey", ""),
            _extract_smiles(props),
        )
    return None


def query_names_batch(name_list):
    """Query PubChem by compound names in batch via POST.

    Returns dict: {name: (cid, inchikey, smiles)} for found names.
    """
    if not name_list:
        return {}
    return _query_names_batch_recursive(name_list)


def _query_names_batch_recursive(name_list):
    """Try batch name lookup; on 404, split and retry recursively."""
    if not name_list:
        return {}

    if len(name_list) == 1:
        name = name_list[0]
        url = f"{PUBCHEM_BASE}/compound/name/property/InChIKey,IsomericSMILES/JSON"
        data = pubchem_request("POST", url, data={"name": name})
        if data and "PropertyTable" in data:
            prop_list = data["PropertyTable"]["Properties"]
            if not prop_list:
                return {}
            props = prop_list[0]
            return {name: (
                str(props.get("CID", "")),
                props.get("InChIKey", ""),
                _extract_smiles(props),
            )}
        return {}

    url = f"{PUBCHEM_BASE}/compound/name/property/InChIKey,IsomericSMILES/JSON"
    joined = "\n".join(name_list)
    data = pubchem_request("POST", url, data={"name": joined})

    if data and "PropertyTable" in data:
        props_list = data["PropertyTable"]["Properties"]
        results = {}
        for i, props in enumerate(props_list):
            if i < len(name_list):
                results[name_list[i]] = (
                    str(props.get("CID", "")),
                    props.get("InChIKey", ""),
                    _extract_smiles(props),
                )
        return results

    # Batch failed (likely 404 due to unknown name) — split in half
    mid = len(name_list) // 2
    left = _query_names_batch_recursive(name_list[:mid])
    right = _query_names_batch_recursive(name_list[mid:])
    left.update(right)
    return left


def query_cid_properties(cid):
    """Fetch InChIKey, IsomericSMILES, InChI for a PubChem CID."""
    url = (f"{PUBCHEM_BASE}/compound/cid/{cid}"
           f"/property/InChIKey,IsomericSMILES,InChI/JSON")
    data = pubchem_request("GET", url)
    if data and "PropertyTable" in data:
        prop_list = data["PropertyTable"]["Properties"]
        if not prop_list:
            return None
        props = prop_list[0]
        return {
            "smiles": _extract_smiles(props),
            "inchi": props.get("InChI", ""),
            "inchikey": props.get("InChIKey", ""),
        }
    return None


CID_BATCH_SIZE = 100


def query_cid_properties_batch(cid_list):
    """Fetch properties for multiple CIDs in a single request.

    PubChem supports comma-separated CIDs in the URL. Returns dict:
        {cid_str: {smiles, inchi, inchikey}} for each found CID.
    """
    if not cid_list:
        return {}
    cid_str = ",".join(str(c) for c in cid_list)
    url = (f"{PUBCHEM_BASE}/compound/cid/{cid_str}"
           f"/property/InChIKey,IsomericSMILES,InChI/JSON")
    data = pubchem_request("GET", url)
    results = {}
    if data and "PropertyTable" in data:
        for props in data["PropertyTable"]["Properties"]:
            cid = str(props.get("CID", ""))
            if cid:
                results[cid] = {
                    "smiles": _extract_smiles(props),
                    "inchi": props.get("InChI", ""),
                    "inchikey": props.get("InChIKey", ""),
                }
    return results
