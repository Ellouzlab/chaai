from dataclasses import dataclass, field
from chaai._core import genecall as _cs


@dataclass
class ORF:
    contig: str
    frame: int
    nuc_start: int
    aa_len: int
    protein: str

    @property
    def strand(self):
        return "+" if self.frame < 3 else "-"


@dataclass
class GenomeDb:
    name: str
    nuc_length: int
    contigs: list = field(default_factory=list)
    orfs: list = field(default_factory=list)

    @property
    def proteins(self):
        return [o.protein for o in self.orfs]

    def write(self, path):
        segs = [(o.frame, o.nuc_start, o.aa_len, o.protein) for o in self.orfs]
        _cs.write_chaai(path, self.name, self.nuc_length,
                        [(c, o, l) for c, o, l in self.contigs], segs)
        return path

    def to_gff(self, path=None):
        lines = ["##gff-version 3"]
        off2contig = sorted(self.contigs, key=lambda c: c[1])
        for i, o in enumerate(self.orfs, 1):
            cname, coff = contig_of(off2contig, o.nuc_start)
            start = o.nuc_start - coff + 1
            end = start + 3 * o.aa_len - 1
            lines.append("\t".join([cname, "chaai", "CDS", str(start), str(end),
                                    ".", o.strand, "0", f"ID=orf_{i}"]))
        text = "\n".join(lines) + "\n"
        if path is not None:
            with open(path, "w") as fh:
                fh.write(text)
        return text

    def to_faa(self, path=None):
        off2contig = sorted(self.contigs, key=lambda c: c[1])
        cidx = {c[0]: k for k, c in enumerate(off2contig, 1)}
        count = {}
        lines = []
        for o in sorted(self.orfs, key=lambda o: o.nuc_start):
            cname, coff = contig_of(off2contig, o.nuc_start)
            k = count[cname] = count.get(cname, 0) + 1
            start = o.nuc_start - coff + 1
            end = start + 3 * o.aa_len - 1
            strand = 1 if o.frame < 3 else -1
            lines.append(f">{cname}_{k} # {start} # {end} # {strand} # ID={cidx[cname]}_{k}")
            lines.append(o.protein)
        text = "\n".join(lines) + "\n"
        if path is not None:
            with open(path, "w") as fh:
                fh.write(text)
        return text

    @classmethod
    def read(cls, path):
        from chaai.db.format import read_chaai
        return read_chaai(path)


def contig_of(sorted_contigs, pos):
    hit = sorted_contigs[0]
    for c in sorted_contigs:
        if c[1] <= pos:
            hit = c
        else:
            break
    return hit[0], hit[1]
