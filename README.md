# chaai

chaai computes average amino-acid identity (AAI) between genomes, from compact viral genomes to large
eukaryotic assemblies. It works in two stages:

1. **`build_db`** translates each genome once into a small database of filtered open reading frames
   (ORFs). A fast built-in ORF filter keeps likely coding ORFs from a six-frame translation;
   Pyrodigal can be used instead.
2. **`dist`** compares every query genome with every reference by seed-chain-align: sparse spaced
   amino-acid seeds are matched, chained within ORF pairs, reciprocal best hits are chosen from the
   chains, and only those pairs are aligned. It reports AAI, genome coverage and a composite score.
   **`search`** avoids the full matrix: a marker index picks a few candidate references per query,
   and only those are compared.

## Install

Requires Python 3.10 or newer and a C++17 compiler with OpenMP.

```bash
pip install .
pip install ".[pyrodigal]"          # optional, for --method pyrodigal
```

Or with conda/mamba, which also supplies the compiler and Pyrodigal:

```bash
mamba env create -f environment.yml
mamba activate chaai
```

The extensions are compiled with `-march=native`, so they run on the machine that built them.
Set `CHAAI_NATIVE=0` before installing for a portable build. On macOS the system clang has no
OpenMP; use the conda route or `brew install libomp`.

## Quick start

```bash
chaai build_db --infolder genomes/ --outfolder db/

# every pair in one set (the upper triangle), or every ordered pair
chaai dist -q db/ -r db/ --small-normal --output dist.tsv
chaai dist -q db/ -r db/ --small-normal --no-triangle --output dist.tsv

# queries against references
chaai dist -q queries_db/ -r refs_db/ --med-normal --output dist.tsv

# the 10 closest references of every query, without the full matrix
chaai search -q queries_db/ -r refs_db/ --med-normal -n 10 --output hits.tsv
```

Pick the preset for your genomes (see [Presets](#presets)). Every command writes a log to
`./chaai_logs`. `dist` and `search` keep a seed cache in `./chaai_temp` (or `--temp-dir`), so
repeated runs on the same databases skip seed generation.

## Commands

### build_db

| option | default | meaning |
|---|---|---|
| `--infolder` / `--infasta` | | a folder with one genome per FASTA file, or one multi-FASTA with one genome per record |
| `--outfolder` | `called_fastas` | where the `.chaai` databases are written, one per genome |
| `--method` | `chaai` | `chaai` (built-in ORF filter) or `pyrodigal` (metagenomic mode; recommended for prokaryotes and their viruses) |
| `--min-orf` | 30 | minimum ORF length in codons |
| `--window` | 300 | window, in nucleotides, of the non-coding background |
| `--no-coding-filter` | off | skip the ORF filter and keep every six-frame ORF of at least `--min-orf` codons |
| `--gff-out` | off | also write the ORFs as GFF3 |
| `--faa` | off | also write the proteins as FASTA with Prodigal-style headers (`>contig_1 # start # end # strand # ID=1_1`); an export only, not usable as a database |
| `--threads` | all cores | |

### dist

| option | default | meaning |
|---|---|---|
| `-q`, `-r` | | query and reference database folders; the same folder compares a set with itself |
| preset flag | | e.g. `--small-normal`; sets the seed and the filters below that are left at their defaults |
| `--no-triangle` | off | with the same folder for `-q` and `-r`, compare every ordered pair instead of the upper triangle |
| `-c`, `--scaled` | preset | keep 1 in *c* seeds |
| `-g`, `--gene-scaled` | 1 | keep 1 in *g* ORFs per genome; matched ORFs and coverage are rescaled by *g*² |
| `--min-aai` | 0.01 | minimum AAI reported |
| `--min-orf-cov` | 0.2, or the preset's | minimum fraction of an ORF covered by its chain |
| `--min-chain-seeds` | 2, or the preset's | minimum anchors per chain |
| `--max-evalue` | 0.1, or the preset's | maximum E-value of an ORF alignment |
| `--xdrop` | 15 | X-drop threshold |
| `--min-matched-orfs` | 0 | minimum reciprocal-best-hit ORF pairs reported |
| `--aggressive-filter` | off | reject low-complexity alignments |
| `--max-pair-matches` | 32,000,000 | cap on anchors held per thread for one genome pair |
| `--temp-dir` | `./chaai_temp` | seed cache location |
| `--output` | `db_distances.tsv` | |
| `--threads` | all cores | |

Without a preset, `dist` uses the med-normal seed shape at `-c 10` with the defaults above.

### search

Takes the `dist` options except `--no-triangle`, `--min-matched-orfs`, `--aggressive-filter`,
`--max-pair-matches` and `--gene-scaled`, and adds:

| option | default | meaning |
|---|---|---|
| `-n`, `--top-n` | 1 | hits reported per query |
| `-k`, `--candidates` | 20 | references aligned per query |
| `--marker-factor` | preset | index 1 in *F* of the marker seeds |
| `--marker-shape` | preset family | seed shape of the markers |
| `--max-posting` | 400 | skip markers held by more references than this; bounds the cost of the screen |
| `--batch-refs` | 0 (all) | references indexed per batch; bounds memory for very large reference sets |
| `--containment-only` | off | report the candidates and their marker scores without aligning them |
| `--exclude-self` | off | skip a reference with the same name as the query |

Without a preset, `search` uses `--small-sensitive`.

## Presets

Four families by genome type, each with a sensitive, normal and fast tier; viruses also have an
ultra tier. A sparser seed (larger *c*) is faster but finds fewer divergent protein pairs, which
also shifts AAI upward among the pairs that remain.

| preset | for | seed shape | *c* | search markers, 1 in *F* |
|---|---|---|---|---|
| `--small-ultra` | viruses, deepest | `4:0,1,3,4` | 1 | `6:0,2,3,4,5,6`, 2 |
| `--small-sensitive` | viruses, order and above | `4:0,1,3,4` | 2 | `6:0,2,3,4,5,6`, 2 |
| `--small-normal` | viruses, family and genus | `4:0,1,4,5` | 8 | `6:0,2,3,4,5,6`, 2 |
| `--small-fast` | viruses, genus screen | `6:0,2,3,4,5,6` | 12 | `6:0,2,3,4,5,6`, 2 |
| `--med-sensitive` | bacteria and archaea | `5:0,1,3,4,5` | 2 | `8:0,1,2,3,4,5,6,9`, 16 |
| `--med-normal` | bacteria and archaea | `6:0,2,3,4,5,6` | 8 | `8:0,1,2,3,4,5,6,9`, 16 |
| `--med-fast` | bacteria and archaea, screen | `6:0,1,2,4,5,6` | 30 | `8:0,1,2,3,4,5,6,9`, 16 |
| `--large-sensitive` | fungi and protists | `6:0,2,3,4,6,7` | 2 | `10:0,…,9`, 256 |
| `--large-normal` | fungi and protists | `6:0,2,3,4,5,6` | 8 | `10:0,…,9`, 256 |
| `--large-fast` | fungi and protists, screen | `6:0,2,3,4,5,6` | 20 | `10:0,…,9`, 256 |
| `--xl-sensitive` | plants and animals | `6:0,2,3,4,6,7` | 4 | `10:0,…,9`, 1024 |
| `--xl-normal` | plants and animals | `8:0,1,2,3,4,5,6,9` | 8 | `10:0,…,9`, 1024 |
| `--xl-fast` | plants and animals, screen | `7:0,2,3,4,5,6,7` | 20 | `10:0,…,9`, 1024 |

A shape `w:positions` lists the care positions of a spaced seed of weight *w*: `6:0,2,3,4,6,7`
spans eight residues and ignores the second and sixth. Every preset sets an E-value cutoff of
0.001. The minimum ORF chain coverage is 0 for the small presets and med-sensitive and 0.1 for the
rest, and small-ultra requires three anchors per chain.

## Output

`dist`, one row per genome pair that passes the filters:

| column | meaning |
|---|---|
| `query`, `ref` | genome names |
| `aai` | identical residues over the aligned reciprocal-best-hit ORF pairs, divided by the summed length of the shorter ORF of each pair |
| `chaining_genome_coverage_q`, `_r` | fraction of each genome's length covered by alignments; the unsuffixed column is the larger |
| `avg_orf_chain_coverage_q`, `_r` | fraction of the matched ORFs' length covered by alignments, weighted by ORF length; the unsuffixed column is the larger |
| `composite_score` | √(AAI × C), where C is the larger of the two genomes' aligned fractions of their coding span |
| `matched_orfs` | aligned reciprocal-best-hit ORF pairs |
| `orf_hit` | `matched_orfs` over the ORF count, the larger of the two genomes |
| `query_len`, `ref_len` | genome lengths in nucleotides |

`search` adds `shared_markers`, `marker_excess` (shared markers above background) and
`marker_containment` (that excess over the smaller marker set, less the reference's typical value)
after `query` and `ref`, followed by the `dist` columns. With `--containment-only` the rows end with
the marker counts of both genomes instead.

## How it works

### ORF filtering

1. In all six frames, every run of at least 30 codons between stop codons is a candidate ORF.
2. Each ORF gets a three-base-periodicity signal-to-noise ratio (SNR), a positional base
   composition score Ψ, a ribosome binding site score and its length. ORFs with SNR ≥ 2.0, at least
   100 codons, and Ψ ≥ 0.015 or a binding site train a hexamer log-odds table against a background
   of 300 bp windows with SNR < 1.6. The table gives every ORF a fourth score, *d*, and is refit once.
3. ORFs are ranked by the sum of the z-scores of *d*, SNR, Ψ and ln *L*.
4. A coding budget is estimated as the ORF length in excess of what chance produces. The stop-codon
   rate *q* is fitted to the ORFs of 30 to 90 codons; chance predicts 2*Mq*²(1−*q*)<sup>*L*</sup>
   ORFs of *L* codons in a contig of *M* bases, and the surplus over all lengths is the budget. ORFs
   are accepted from the top of the ranking until it is spent. Of two ORFs overlapping by more than
   60% of the shorter, the one with the lower z(Ψ) + z(ln *L*) + z(*d*) is discarded.
5. Each kept ORF is trimmed to its first ATG, GTG or TTG within 30 codons that leaves at least 30
   codons, translated and packed at 5 bits per residue.

The filter does not model splicing and is not a gene annotator. It aims to keep the protein
content of a genome while removing most of the spurious ORFs a six-frame translation produces; the
alignment step tolerates the errors that remain.

### Distance

1. **Seeds.** Spaced amino-acid seeds are taken from the packed ORFs, and one in *c* is kept by
   FracMinHash, so the same seeds survive in every genome. A seed occurring more often than its
   genome's composition predicts, above *E* + 4√*E*, is dropped. Optionally, one in *g* ORFs is kept
   the same way before seeding.
2. **Anchors.** The sorted seed lists of two genomes are merge-joined. A seed found *n*<sub>q</sub>
   and *n*<sub>r</sub> times gives at most 8 × 8 anchors, each weighted
   1/√(*n*<sub>q</sub>*n*<sub>r</sub>). A pair with fewer anchors than one chain needs is dropped.
3. **Chaining.** Within each ORF pair, anchors are chained by dynamic programming: an anchor links
   to the best earlier anchor whose diagonal is within 25 residues, looking back up to 500 residues.
   One chain is taken per predecessor tree. Each ORF pair is scored by its best chain, and only
   reciprocal best hits go on.
4. **Alignment.** A chain is split where its anchor diagonals spread more than 12 residues. Each
   segment is extended by ungapped X-drop (X = 15, BLOSUM62) along the diagonals of its end anchors,
   then aligned by banded Smith–Waterman (BLOSUM62, gap open 11, extend 1) around the median anchor
   diagonal. Local maxima give high-scoring segment pairs, filtered by E-value.
5. **Scores.** AAI is the identical residues over all aligned pairs divided by the summed length of
   the shorter ORF of each pair. Coverage C is the larger of the two genomes' aligned fractions of
   their coding span, so a partial assembly is not penalised, and the composite score is √(C × AAI).
   The composite score discounts a high AAI that rests on a few shared genes.

### Search

Markers are the seeds that pass a second FracMinHash filter, 1 in *F*, with a heavier seed shape
per preset family. Reference markers go into an inverted index. Each query is scored against every
reference it shares markers with, using a background sample of up to 256 evenly spaced references:
once by the shared markers in excess of the query's background rate, and once by that excess over
the smaller marker set, less the reference's typical value. Candidates are taken alternately from
the two rankings, 20 by default, and compared with the full distance pipeline. The reported hits
are the closest references among the candidates, which usually but not always include the
exhaustive best.

## Python API

```python
import chaai

db = chaai.build_db("genome.fna")                     # ORF filter -> GenomeDb
db = chaai.build_db("genome.fna", method="pyrodigal")
db = chaai.build_db(seq="ATG...", name="g")           # from a sequence string

db.orfs                                               # [ORF(contig, frame, nuc_start, aa_len, protein)]
db.proteins                                           # [str]
db.write("genome.chaai")
db.to_gff("genome.gff")
db.to_faa("genome.faa")                               # export only
db2 = chaai.read_db("genome.chaai")

rows = chaai.dist([db, db2], preset="small-normal")   # GenomeDb objects, .chaai paths or a folder
chaai.build_db_folder("genomes/", "db/")

for frame, nuc_start, aa_len, protein in chaai.find_genes("contig", seq):
    ...                                               # frames 0-2 forward, 3-5 reverse
```

## License

GPL-3.0-only. See [LICENSE](LICENSE).

## Contact

Muhammad Sulman, sulmanmu40@gmail.com

Oualid Ellouz, oualid.ellouz@agr.gc.ca
