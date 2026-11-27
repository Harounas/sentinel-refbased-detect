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

    def efetch_fasta(self, accessions, attempts=5, chunk=200):
        """FASTA for a list of accessions, in batches and by POST.

        The id list goes in the request body, not the query string: a few
        hundred accessions make a URL long enough for NCBI to answer 414
        Request-URI Too Long, which it does silently as far as the caller is
        concerned. Batching keeps each request well inside the limit and is
        also what NCBI asks for above a couple of hundred ids.
        """
        accessions = list(accessions)
        out = []
        for i in range(0, len(accessions), chunk):
            batch = accessions[i:i + chunk]
            params = {"db": "nuccore", "id": ",".join(batch),
                      "rettype": "fasta", "retmode": "text",
                      "tool": TOOL, "email": self.email}
            if self.api_key:
                params["api_key"] = self.api_key
            data = urllib.parse.urlencode(params).encode()
            last = None
            for attempt in range(attempts):
                self._throttle()
                try:
                    req = urllib.request.Request(EUTILS + "efetch.fcgi", data=data)
                    with urllib.request.urlopen(req, timeout=300) as resp:
                        text = resp.read().decode()
                    if text.lstrip().startswith(">"):
                        out.append(text if text.endswith("\n") else text + "\n")
                        last = None
                        break
                    last = RuntimeError("no FASTA returned: " + text[:200])
                except Exception as exc:  # noqa: BLE001
                    last = exc
                time.sleep(2 ** attempt)
            if last is not None:
                raise last
        return "".join(out)
