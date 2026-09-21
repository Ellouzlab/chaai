import os
import csv
import glob
import shutil
import tempfile

from chaai._core import dist as _cd
from chaai.cli.presets import DEFAULT_PRESET, PRESETS
from chaai.db.model import GenomeDb


def as_paths(items, tmpdir):
    if isinstance(items, str):
        if os.path.isdir(items):
            return sorted(glob.glob(os.path.join(items, "*.chaai")))
        return [items]
    paths = []
    for i, it in enumerate(items):
        if isinstance(it, GenomeDb):
            p = os.path.join(tmpdir, f"{i}_{it.name}.chaai")
            it.write(p); paths.append(p)
        else:
            paths.append(it)
    return paths


def dist(queries, refs=None, threads=1, is_self=None, min_aai=0.01, preset=DEFAULT_PRESET, **kw):
    cfg = PRESETS[preset]
    kw.setdefault("seed_shape", cfg["shape"])
    kw.setdefault("seed_scaled", cfg["scaled"])
    tmp = tempfile.mkdtemp(prefix="chaai_dist_")
    try:
        q = as_paths(queries, tmp)
        r = q if refs is None else as_paths(refs, tmp)
        if is_self is None:
            is_self = refs is None
        out = os.path.join(tmp, "dist.tsv")
        _cd.run_dist(q, r, out, threads, is_self, min_aai=min_aai, **kw)
        with open(out) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        for row in rows:
            for k, v in row.items():
                if k not in ("query", "ref"):
                    try:
                        row[k] = float(v)
                    except ValueError:
                        pass
        return rows
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def run_dist_files(queries, refs, output_path, threads, is_self, **opts):
    return _cd.run_dist(queries, refs, output_path, threads, is_self, **opts)
