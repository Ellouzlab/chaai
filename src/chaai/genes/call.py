import os
import glob
import tempfile
import shutil

from chaai._core import genecall as _cs
from chaai.db.alphabet import sanitize_protein
from chaai.db.fasta import read_fasta
from chaai.db.format import read_chaai
from chaai.db.model import ORF, GenomeDb


def gene_finder_params(min_orf, coding_filter, map_window, threads):
    p = _cs.GeneFinderParams()
    p.min_orf, p.coding_filter, p.map_window, p.threads = \
        min_orf, coding_filter, map_window, threads
    return p


def call_chaai(contigs, params):
    if len(contigs) == 1:
        cname, seq = contigs[0]
        return [ORF(cname, fr, ns, al, pr) for fr, ns, al, pr in _cs.find_genes(cname, seq, params)]
    tmp = tempfile.mkdtemp(prefix="chaai_build_")
    try:
        with open(os.path.join(tmp, "g.fna"), "w") as fh:
            for c, s in contigs:
                fh.write(">" + c + "\n")
                for i in range(0, len(s), 80):
                    fh.write(s[i:i + 80] + "\n")
        out = os.path.join(tmp, "out"); os.makedirs(out)
        _cs.build_from_folder(tmp, out, params.threads, params.coding_filter,
                               params.min_orf, params.map_window, False)
        return read_chaai(glob.glob(os.path.join(out, "*.chaai"))[0]).orfs
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def call_pyrodigal(contigs, meta=True, min_orf=None):
    import pyrodigal
    gf = pyrodigal.GeneFinder(meta=meta)
    orfs = []
    offset = 0
    for cname, seq in contigs:
        if not meta and len(seq) >= 20000:
            gf = pyrodigal.GeneFinder(); gf.train(seq)
        genes = gf.find_genes(seq.encode()) if isinstance(seq, str) else gf.find_genes(seq)
        for g in genes:
            prot = sanitize_protein(g.translate())
            aa_len = len(prot)
            if min_orf and aa_len < min_orf:
                continue
            frame = (g.begin - 1) % 3 if g.strand == 1 else 3 + ((len(seq) - g.begin) % 3)
            nuc_start = offset + (g.begin - 1)
            orfs.append(ORF(cname, frame, nuc_start, aa_len, prot))
        offset += len(seq)
    return orfs


def build_db(source=None, method="chaai", name=None, seq=None,
             min_orf=30, coding_filter=True, map_window=300,
             threads=1):
    if seq is not None:
        contigs = [(name or "genome", seq.upper())]
    elif source is not None:
        contigs = [(c, s.upper()) for c, s in read_fasta(source)]
        name = name or os.path.splitext(os.path.basename(source))[0]
    else:
        raise ValueError("provide source=<fasta path> or seq=<sequence>")
    name = name or "genome"
    nuc_length = sum(len(s) for _, s in contigs)
    contig_table, off = [], 0
    for c, s in contigs:
        contig_table.append((c, off, len(s))); off += len(s)

    if method == "pyrodigal":
        orfs = call_pyrodigal(contigs, min_orf=min_orf)
    elif method == "chaai":
        orfs = call_chaai(contigs, gene_finder_params(min_orf, coding_filter, map_window, threads))
    else:
        raise ValueError(f"unknown method: {method!r}")
    return GenomeDb(name, nuc_length, contig_table, orfs)


def build_db_folder(infolder, outfolder, method="chaai", exts=(".fa", ".fna", ".fasta", ".fas"),
                    gff_out=False, faa_out=False, **kw):
    os.makedirs(outfolder, exist_ok=True)
    out = []
    for path in sorted(glob.glob(os.path.join(infolder, "*"))):
        if os.path.splitext(path)[1].lower() not in exts:
            continue
        stem = os.path.splitext(os.path.basename(path))[0]
        sk = build_db(path, method=method, **kw)
        sk.write(os.path.join(outfolder, stem + ".chaai"))
        if gff_out:
            sk.to_gff(os.path.join(outfolder, stem + ".gff"))
        if faa_out:
            sk.to_faa(os.path.join(outfolder, stem + ".faa"))
        out.append(sk)
    return out
