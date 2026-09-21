import argparse

from chaai import __version__
import sys
from multiprocessing import cpu_count


def verify(arguments):
    if arguments.command is None:
        print("Please specify a command: build_db, dist or search")
        sys.exit(1)

    if arguments.command == "build_db":
        if arguments.infasta is None and arguments.infolder is None:
            print("Please specify either --infasta or --infolder.")
            sys.exit(1)
        if arguments.infasta is not None and arguments.infolder is not None:
            print("Please specify only one of --infasta or --infolder.")
            sys.exit(1)


def _add_seed_args(parser):
    parser.add_argument("-c", "--scaled", type=int, default=10,
                        help="Seed sampling rate: keep 1/c of seeds (default: 10)")


def _add_preset_args(parser):
    from chaai.cli.presets import PRESETS, HELP
    preset = parser.add_mutually_exclusive_group()
    for name, cfg in PRESETS.items():
        preset.add_argument(f"--{name}", dest="preset", action="store_const", const=name,
                            help=HELP.get(name, f"w{cfg['w']} c={cfg['scaled']}"))


def argparser():
    args = argparse.ArgumentParser(
        prog="chaai",
        description="chaai - Chaining Amino Acid k-mers for AAI: gene calling and genome distance"
    )
    args.add_argument('-v', '--version', action='version', version=__version__)
    subparsers = args.add_subparsers(dest='command', help='sub-command help')

    build_parser = subparsers.add_parser('build_db', help='Call genes and write a .chaai database per genome')
    build_parser.add_argument("--infasta", required=False, help="Input FASTA file (each record = one genome)")
    build_parser.add_argument("--infolder", required=False, help="Input folder (each FASTA file = one genome)")
    build_parser.add_argument("--outfolder", default="called_fastas", help="Output folder (default: called_fastas)")
    build_parser.add_argument("--method", choices=["chaai", "pyrodigal"], default="chaai",
                              help="Gene finder (default: chaai)")
    build_parser.add_argument("--gff-out", action="store_true", default=False,
                              help="Also write a GFF of the ORFs next to each database")
    build_parser.add_argument("--faa", action="store_true", default=False,
                              help="Also write the called proteins as a FASTA (.faa) next to each database, with Prodigal-style headers; an export only, it cannot be used as a database")
    build_parser.add_argument("--no-coding-filter", action="store_true", default=False,
                              help="Disable the coding filter (enabled by default)")
    build_parser.add_argument("--window", type=int, default=300,
                              help="Coding map window in nucleotides (default: 300)")
    build_parser.add_argument("--min-orf", type=int, default=30,
                              help="Minimum ORF length in amino acids (default: 30)")
    build_parser.add_argument("--verbose", action="store_true", default=False,
                              help="Print per-genome stats")
    build_parser.add_argument("--threads", type=int, default=cpu_count(), help="Number of threads")

    dist_parser = subparsers.add_parser('dist', help='Calculate distances between .chaai databases')
    dist_parser.add_argument("-q", "--qdbfolder", required=True, help="Query .chaai folder")
    dist_parser.add_argument("-r", "--rdbfolder", required=True, help="Reference .chaai folder")
    _add_seed_args(dist_parser)
    _add_preset_args(dist_parser)
    dist_parser.add_argument("--min-matched-orfs", type=int, default=0,
                             help="Minimum RBH ORF pairs to report (default: 0)")
    dist_parser.add_argument("--min-aai", type=float, default=0.01,
                             help="Minimum AAI to report (default: 0.01)")
    dist_parser.add_argument("--min-orf-cov", type=float, default=0.20,
                             help="Minimum chain coverage of ORF (default: 0.20)")
    dist_parser.add_argument("--min-chain-seeds", type=int, default=2,
                             help="Minimum anchors per chain (default: 2)")
    dist_parser.add_argument("--aggressive-filter", action="store_true", default=False,
                             help="Reject low-complexity chain alignments (repeat filter)")
    dist_parser.add_argument("--xdrop", type=int, default=15,
                             help="X-drop threshold for ungapped extension (default: 15)")
    dist_parser.add_argument("--max-evalue", type=float, default=0.1,
                             help="Maximum E-value for ORF alignments (default: 0.1)")
    dist_parser.add_argument("--temp-dir", default="",
                             help="Where to keep the seed cache (default: chaai_temp here)")
    dist_parser.add_argument("--output", default="db_distances.tsv", help="Output TSV")
    dist_parser.add_argument("--no-triangle", action="store_true", default=False,
                             help="Same query and reference folder: compare every ordered pair instead of the upper triangle")
    dist_parser.add_argument("--max-pair-matches", type=int, default=None,
                             help="Cap on seed matches held per thread for one genome pair, "
                                  "16 B each (default: 32000000; 0=unlimited)")
    dist_parser.add_argument("-g", "--gene-scaled", type=int, default=1,
                             help="Gene sampling rate: keep 1/g of the ORFs of every genome "
                                  "(FracMinHash), matched ORFs and coverage rescaled by g^2 "
                                  "(default: 1 = all genes)")
    dist_parser.add_argument("--threads", type=int, default=cpu_count(), help="Number of threads")

    search_parser = subparsers.add_parser(
        'search', help='Find the closest references: marker index, then dist on the best candidates')
    search_parser.add_argument("-q", "--qdbfolder", required=True, help="Query .chaai folder")
    search_parser.add_argument("-r", "--rdbfolder", required=True, help="Reference .chaai folder")
    search_parser.add_argument("-n", "--top-n", type=int, default=1,
                               help="Hits reported per query (default: 1)")
    search_parser.add_argument("-k", "--candidates", type=int, default=None,
                               help="References aligned per query, chosen by the index (default: per preset)")
    search_parser.add_argument("--marker-factor", type=int, default=None,
                               help="Index 1 in F of the marker seeds (default: per preset)")
    search_parser.add_argument("--marker-shape", default=None,
                               help="Seed shape of the index markers, at the preset's rate (default: per preset family)")
    search_parser.add_argument("--max-posting", type=int, default=400,
                               help="Skip markers held by more than this many references of a batch; bounds screen cost (default: 400, 0 = no cap)")
    search_parser.add_argument("--batch-refs", type=int, default=0,
                               help="References indexed per batch; bounds memory (default: 0, all at once)")
    search_parser.add_argument("--containment-only", action="store_true", default=False,
                               help="Report the index's candidates without aligning them")
    search_parser.add_argument("--exclude-self", action="store_true", default=False,
                               help="Skip a reference with the same name as the query")
    _add_seed_args(search_parser)
    _add_preset_args(search_parser)
    search_parser.add_argument("--min-aai", type=float, default=0.01,
                               help="Minimum AAI to report (default: 0.01)")
    search_parser.add_argument("--min-orf-cov", type=float, default=0.20,
                               help="Minimum chain coverage of ORF (default: 0.20)")
    search_parser.add_argument("--min-chain-seeds", type=int, default=2,
                               help="Minimum anchors per chain (default: 2)")
    search_parser.add_argument("--max-evalue", type=float, default=0.1,
                               help="Maximum E-value for ORF alignments (default: 0.1)")
    search_parser.add_argument("--xdrop", type=int, default=15,
                               help="X-drop threshold for ungapped extension (default: 15)")
    search_parser.add_argument("--temp-dir", default="",
                               help="Where to keep the seed cache (default: chaai_temp here)")
    search_parser.add_argument("--output", default="search_results.tsv", help="Output TSV")
    search_parser.add_argument("--threads", type=int, default=cpu_count())

    arguments = args.parse_args()
    verify(arguments)
    return arguments
