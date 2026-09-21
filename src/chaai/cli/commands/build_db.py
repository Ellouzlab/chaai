import os
import sys
import logging
import time


class BuildDbHandler:

    def __init__(self, args):
        logging.info("=== chaai build_db ===")
        from chaai import build_db, build_db_folder
        from chaai._core.genecall import build_from_fasta, build_from_folder

        out = args.outfolder
        os.makedirs(out, exist_ok=True)
        method = getattr(args, "method", "chaai")
        coding_filter = not getattr(args, "no_coding_filter", False)
        min_orf = getattr(args, "min_orf", 30)
        window = getattr(args, "window", 300)
        gff_out = getattr(args, "gff_out", False)
        faa_out = getattr(args, "faa", False)
        threads = args.threads
        t0 = time.time()

        fast = (method == "chaai" and not gff_out and not faa_out)

        if args.infasta:
            if not os.path.isfile(args.infasta):
                logging.error(f"Input file not found: {args.infasta}"); sys.exit(1)
            if fast:
                n = build_from_fasta(args.infasta, out, threads, coding_filter,
                                      min_orf, window, args.verbose)
            else:
                from chaai.db.fasta import read_fasta
                n = 0
                for cname, seq in read_fasta(args.infasta):
                    sk = build_db(seq=seq, name=cname, method=method, min_orf=min_orf,
                                  coding_filter=coding_filter, map_window=window,
                                  threads=threads)
                    sk.write(os.path.join(out, cname + ".chaai"))
                    if gff_out:
                        sk.to_gff(os.path.join(out, cname + ".gff"))
                    if faa_out:
                        sk.to_faa(os.path.join(out, cname + ".faa"))
                    n += 1
        else:
            if not os.path.isdir(args.infolder):
                logging.error(f"Input folder not found: {args.infolder}"); sys.exit(1)
            if fast:
                n = build_from_folder(args.infolder, out, threads, coding_filter,
                                       min_orf, window, args.verbose)
            else:
                sks = build_db_folder(args.infolder, out, method=method,
                                      gff_out=gff_out, faa_out=faa_out, min_orf=min_orf, coding_filter=coding_filter,
                                      map_window=window, threads=threads)
                n = len(sks)

        logging.info(f"Built {n} databases ({method}) in {time.time() - t0:.1f}s")
        logging.info("=== chaai build_db complete ===")
