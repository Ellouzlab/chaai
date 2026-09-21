import sys
import logging
import time

from chaai.cli.presets import (PRESETS, SEARCH_INDEX, DEFAULT_SEARCH_PRESET,
                               apply_preset, dist_options, list_dbs)


class SearchHandler:
    def __init__(self, args):
        q_paths = list_dbs(args.qdbfolder)
        r_paths = list_dbs(args.rdbfolder)
        logging.info(f"Query: {len(q_paths)}, Ref: {len(r_paths)}")

        try:
            from chaai._core.search import run_search
        except ImportError as e:
            logging.error(f"Cannot import search module: {e}")
            sys.exit(1)

        if not getattr(args, 'preset', None):
            args.preset = DEFAULT_SEARCH_PRESET
        apply_preset(args)
        opts = dist_options(args)

        factor, k, mshape = SEARCH_INDEX[args.preset]
        if args.marker_factor is not None:
            factor = args.marker_factor
        if args.candidates is not None:
            k = args.candidates
        if args.marker_shape:
            mshape = args.marker_shape
        k = max(k, args.top_n)
        logging.info(f"Preset {args.preset}: seeds {opts['seed_shape']} c={opts['seed_scaled']}, "
                     f"index markers {mshape or opts['seed_shape']} thinned 1/{factor}, "
                     f"{k} candidates aligned per query")

        t0 = time.time()
        n = run_search(
            q_paths, r_paths, args.output, args.threads,
            top_n=args.top_n, candidates=k,
            seed_shape=opts['seed_shape'], seed_scaled=opts['seed_scaled'],
            marker_shape=mshape or "", marker_factor=factor, max_posting=args.max_posting,
            batch_refs=args.batch_refs,
            containment_only=args.containment_only, exclude_self=args.exclude_self,
            max_seed_occ=opts['max_seed_occ'], chain_band=opts['chain_band'],
            chain_max_qgap=opts['chain_max_qgap'], min_chain_seeds=opts['min_chain_seeds'],
            min_aai=opts['min_aai'], min_matched_orfs=opts['min_matched_orfs'],
            min_orf_cov=opts['min_orf_cov'], aggressive_filter=opts['aggressive_filter'],
            xdrop=opts['xdrop'], max_evalue=opts['max_evalue'],
            max_pair_matches=opts['max_pair_matches'], gene_scaled=opts['gene_scaled'],
            temp_dir=opts['temp_dir'],
        )
        dt = time.time() - t0
        logging.info(f"{n} hits in {dt:.1f}s")
        logging.info("=== done ===")
