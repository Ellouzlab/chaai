import struct

from chaai.db.alphabet import unpack_5bit
from chaai.db.model import ORF, GenomeDb, contig_of


def read_chaai(path):
    data = open(path, "rb").read()
    off = 0
    magic = data[off:off + 8]; off += 8
    if magic != b"CHAAIDB2":
        if magic == b"CHAAIDB\x00":
            raise ValueError(f"built by an older chaai with 32-bit coordinates, rebuild: {path}")
        raise ValueError(f"not a .chaai file: {path}")
    (nl,) = struct.unpack_from("<H", data, off); off += 2
    name = data[off:off + nl].decode(); off += nl
    (nuc_length,) = struct.unpack_from("<Q", data, off); off += 8
    off += 8
    contigs = []
    (nc,) = struct.unpack_from("<I", data, off); off += 4
    for _ in range(nc):
        (cl,) = struct.unpack_from("<H", data, off); off += 2
        cn = data[off:off + cl].decode(); off += cl
        a, b = struct.unpack_from("<QQ", data, off); off += 16
        contigs.append((cn, a, b))
    (ns,) = struct.unpack_from("<Q", data, off); off += 8
    segs = []
    for _ in range(ns):
        frame = data[off]; off += 1
        (nuc_start,) = struct.unpack_from("<Q", data, off); off += 8
        (aa_len,) = struct.unpack_from("<H", data, off); off += 2
        segs.append((frame, nuc_start, aa_len))
    scontigs = sorted(contigs, key=lambda c: c[1]) or [(name, 0, nuc_length)]
    orfs = []
    p = off
    for frame, nuc_start, aa_len in segs:
        prot = unpack_5bit(data, p, aa_len)
        p += (aa_len * 5 + 7) // 8
        cname, _ = contig_of(scontigs, nuc_start)
        orfs.append(ORF(cname, frame, nuc_start, aa_len, prot))
    return GenomeDb(name, nuc_length, contigs, orfs)


def read_db(path):
    return GenomeDb.read(path)
