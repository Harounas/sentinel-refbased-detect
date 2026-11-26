"""Task fetch_references, with a run-wide cache so a reference shared by many
samples is downloaded from NCBI once."""

import csv
import hashlib
import shutil
from pathlib import Path

from .utils import PipelineError, log


def ref_tag(row):
    tag = "{}__{}__{}".format(row["taxid"], row["virus_code"], row["ref_index"])
    return "".join(c if (c.isalnum() or c in "._-") else "_" for c in tag)


def fetch_references(candidates_tsv, refs_dir, cache_dir, entrez):
    """Write refs_dir/<taxid>__<code>__<index>.fasta for every candidate row;
    return the list of paths."""
    refs_dir, cache_dir = Path(refs_dir), Path(cache_dir)
    refs_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    written = []
    with open(candidates_tsv) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            accs = [a for a in row["accessions"].split(":") if a]
            tag = ref_tag(row)
            path = refs_dir / (tag + ".fasta")
            key = hashlib.sha1(":".join(accs).encode()).hexdigest()[:16]
            cached = cache_dir / (key + ".fasta")
            if not cached.exists():
                text = entrez.efetch_fasta(accs)
                n_seq = text.count(">")
                if n_seq != len(accs):
                    raise PipelineError("%s: asked for %d sequences, NCBI returned %d"
                                        % (tag, len(accs), n_seq))
                tmp = cached.with_name(cached.name + ".part")
                tmp.write_text(text if text.endswith("\n") else text + "\n")
                tmp.rename(cached)
            shutil.copyfile(cached, path)
            written.append(path)
    log.debug("wrote %d reference FASTA files to %s", len(written), refs_dir)
    return written
