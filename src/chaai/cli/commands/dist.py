import os
import logging
import time

from chaai.aai import run_dist_files
from chaai.cli.presets import apply_preset, dist_options, list_dbs


class DistHandler:
    def __init__(self, args):
        q_paths = list_dbs(args.qdbfolder)
        r_paths = list_dbs(args.rdbfolder)

        logging.info(f"Query: {len(q_paths)}, Ref: {len(r_paths)}")
        is_self = os.path.realpath(args.qdbfolder) == os.path.realpath(args.rdbfolder)
        if args.no_triangle and not is_self:
            raise SystemExit("--no-triangle needs the same query and reference folder")
        if is_self:
            logging.info("Self-mode: every ordered pair" if args.no_triangle else "Self-mode: upper triangle only")

        apply_preset(args)

        t0 = time.time()
        n = run_dist_files(q_paths, r_paths, args.output, args.threads, is_self,
                           no_triangle=args.no_triangle, **dist_options(args))
        dt = time.time() - t0
        logging.info(f"{n} pairs in {dt:.1f}s" + (f" ({n/dt:.0f}/s)" if n and dt else ""))
        logging.info("=== done ===")
