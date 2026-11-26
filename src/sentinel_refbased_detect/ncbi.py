"""NCBI E-utilities access with retries and rate limiting."""

import json
import threading
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOOL = "sentinel_refbased_detect"

_lock = threading.Lock()
_last = [0.0]


class Entrez:
    def __init__(self, email, api_key=None):
        if not email:
            raise ValueError("an NCBI e-mail is required (--ncbi-email or $NCBI_EMAIL)")
        self.email = email
        self.api_key = api_key or ""
        # 10 requests/s with a key, 3/s without.
        self.gap = 0.12 if self.api_key else 0.4

    def _throttle(self):
        with _lock:
            wait = _last[0] + self.gap - time.time()
            if wait > 0:
                time.sleep(wait)
            _last[0] = time.time()

    def _url(self, endpoint, params):
        params = dict(params, tool=TOOL, email=self.email)
        if self.api_key:
            params["api_key"] = self.api_key
        return EUTILS + endpoint + "?" + urllib.parse.urlencode(params)

    def json(self, endpoint, attempts=4, **params):
        url = self._url(endpoint, dict(params, retmode="json"))
        last = None
        for attempt in range(attempts):
            self._throttle()
            try:
                with urllib.request.urlopen(url, timeout=120) as fh:
                    return json.load(fh)
            except Exception as exc:  # noqa: BLE001
                last = exc
                time.sleep(2 ** attempt)
        raise last

    def efetch_fasta(self, accessions, attempts=5):
        url = self._url("efetch.fcgi", {"db": "nuccore", "id": ",".join(accessions),
                                        "rettype": "fasta", "retmode": "text"})
        last = None
        for attempt in range(attempts):
            self._throttle()
            try:
                with urllib.request.urlopen(url, timeout=120) as resp:
                    text = resp.read().decode()
                if text.lstrip().startswith(">"):
                    return text
                last = RuntimeError("no FASTA returned: " + text[:200])
            except Exception as exc:  # noqa: BLE001
                last = exc
            time.sleep(2 ** attempt)
        raise last
