"""Small helpers: subprocess pipelines, tool checks, file localization."""

import gzip
import logging
import os
import shutil
import subprocess
import time
import urllib.request
from pathlib import Path

log = logging.getLogger("sentinel_refbased_detect")

CORE_TOOLS = ["kraken2", "bwa", "samtools", "ivar"]


class PipelineError(RuntimeError):
    pass


def _fmt(cmd):
    return " ".join(str(c) for c in cmd)


def run(cmd, *, stdout=None, stderr=None, cwd=None, check=True):
    """Run one command; stdout/stderr may be paths or open handles."""
    log.debug("$ %s", _fmt(cmd))
    out_fh = open(stdout, "w") if isinstance(stdout, (str, Path)) else stdout
    err_fh = open(stderr, "a") if isinstance(stderr, (str, Path)) else stderr
    try:
        res = subprocess.run([str(c) for c in cmd], stdout=out_fh,
                             stderr=err_fh, cwd=cwd)
    finally:
        if out_fh is not stdout:
            out_fh.close()
        if err_fh is not stderr:
            err_fh.close()
    if check and res.returncode != 0:
        raise PipelineError("command failed (rc=%d): %s" % (res.returncode, _fmt(cmd)))
    return res.returncode


def pipeline(cmds, *, stdout=None, stderr=None, cwd=None, check=True):
    """Run cmd1 | cmd2 | ... and return the return code of the LAST command.

    With check=True every stage must succeed (like `set -o pipefail`).
    """
    log.debug("$ %s", " | ".join(_fmt(c) for c in cmds))
    out_fh = open(stdout, "w") if isinstance(stdout, (str, Path)) else stdout
    err_fh = open(stderr, "a") if isinstance(stderr, (str, Path)) else stderr
    procs = []
    try:
        prev = None
        for i, cmd in enumerate(cmds):
            last = i == len(cmds) - 1
            p = subprocess.Popen([str(c) for c in cmd], stdin=prev,
                                 stdout=out_fh if last else subprocess.PIPE,
                                 stderr=err_fh, cwd=cwd)
            if prev is not None:
                prev.close()
            prev = p.stdout
            procs.append(p)
        rcs = [p.wait() for p in procs]
    finally:
        if out_fh is not stdout and out_fh is not None:
            out_fh.close()
        if err_fh is not stderr and err_fh is not None:
            err_fh.close()
    if check and any(rcs):
        raise PipelineError("pipeline failed (rc=%s): %s"
                            % (rcs, " | ".join(_fmt(c) for c in cmds)))
    return rcs[-1]


def capture(cmd, cwd=None):
    log.debug("$ %s", _fmt(cmd))
    res = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, cwd=cwd)
    if res.returncode != 0:
        raise PipelineError("command failed (rc=%d): %s\n%s"
                            % (res.returncode, _fmt(cmd), res.stderr[-2000:]))
    return res.stdout


def missing_tools(names):
    return [n for n in names if shutil.which(n) is None]


def require_tools(names):
    miss = missing_tools(names)
    if miss:
        raise PipelineError(
            "required tool(s) not on PATH: %s. Install them, e.g. "
            "`mamba env create -f environment.yml`." % ", ".join(miss))


def localize(src, dest_dir, name=None):
    """Return a local path for src: a local file, gs:// URI or http(s) URL."""
    src = str(src)
    if not (src.startswith("gs://") or src.startswith("http://")
            or src.startswith("https://")):
        p = Path(src).expanduser()
        if not p.exists():
            raise PipelineError("file not found: %s" % src)
        return p.resolve()
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / (name or src.rstrip("/").split("/")[-1])
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_name(dest.name + ".part")
    if src.startswith("gs://"):
        tool = shutil.which("gcloud")
        if tool:
            run([tool, "storage", "cp", src, tmp])
        elif shutil.which("gsutil"):
            run(["gsutil", "cp", src, tmp])
        else:
            raise PipelineError("need gcloud or gsutil to fetch %s" % src)
    else:
        download(src, tmp)
    tmp.rename(dest)
    return dest


def download(url, dest, attempts=4):
    last = None
    for attempt in range(attempts):
        try:
            log.info("downloading %s", url)
            urllib.request.urlretrieve(url, str(dest))
            return
        except Exception as exc:  # noqa: BLE001 - retry anything
            last = exc
            time.sleep(2 ** attempt)
    raise PipelineError("download failed: %s: %s" % (url, last))


def open_text(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def write_done(directory):
    Path(directory, ".done").write_text(time.strftime("%Y-%m-%dT%H:%M:%S\n"))


def is_done(directory):
    return Path(directory, ".done").exists()


def file_nonempty(path):
    return os.path.exists(path) and os.path.getsize(path) > 0
