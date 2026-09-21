import os
import sys
import logging
from typing import List, Tuple


def _p(w, shape, scaled, min_orf_cov, budget=32_000_000, min_chain_seeds=2):
    return {'w': w, 'shape': shape, 'scaled': scaled, 'min_orf_cov': min_orf_cov,
            'max_pair_matches': budget,
            'min_aai': 0.01, 'max_evalue': 0.001, 'min_chain_seeds': min_chain_seeds,
            'xdrop': 15}


_S6='6:0,2,3,4,5,6'
_S6F='6:0,1,2,4,5,6'
_S7='7:0,2,3,4,5,6,7'
_S8='8:0,1,2,3,4,5,6,9'

PRESETS = {
    'small-ultra':     _p(4, '4:0,1,3,4', 1, 0,   32_000_000, min_chain_seeds=3),
    'small-sensitive': _p(4, '4:0,1,3,4', 2, 0,   32_000_000),
    'small-normal':    _p(4, '4:0,1,4,5', 8, 0,    8_000_000),
    'small-fast':      _p(6, '6:0,2,3,4,5,6', 12, 0, 2_000_000),
    'med-sensitive':   _p(5, '5:0,1,3,4,5', 2, 0,  32_000_000),
    'med-normal':      _p(6, _S6,         8, 0.1,  8_000_000),
    'med-fast':        _p(6, _S6F,       30, 0.1,  2_000_000),
    'large-sensitive': _p(6, '6:0,2,3,4,6,7', 2, 0.1, 32_000_000),
    'large-normal':    _p(6, _S6,         8, 0.1,  8_000_000),
    'large-fast':      _p(6, _S6,        20, 0.1,  2_000_000),
    'xl-sensitive':    _p(6, '6:0,2,3,4,6,7', 4, 0.1, 32_000_000),
    'xl-normal':       _p(8, _S8,         8, 0.1,  8_000_000),
    'xl-fast':         _p(7, _S7,        20, 0.1,  2_000_000),
}
DEFAULT_SEARCH_PRESET = 'small-sensitive'
_MK = {'small': '6:0,2,3,4,5,6', 'med': '8:0,1,2,3,4,5,6,9', 'large': '10:0,1,2,3,4,5,6,7,8,9', 'xl': '10:0,1,2,3,4,5,6,7,8,9'}
SEARCH_INDEX = {
    'small-ultra':     (2, 20, _MK['small']),
    'small-sensitive': (2, 20, _MK['small']),
    'small-normal':    (2, 20, _MK['small']),
    'small-fast':      (2, 20, _MK['small']),
    'med-sensitive':   (16, 20, _MK['med']),
    'med-normal':      (16, 20, _MK['med']),
    'med-fast':        (16, 20, _MK['med']),
    'large-sensitive': (256, 20, _MK['large']),
    'large-normal':    (256, 20, _MK['large']),
    'large-fast':      (256, 20, _MK['large']),
    'xl-sensitive':    (1024, 20, _MK['xl']),
    'xl-normal':       (1024, 20, _MK['xl']),
    'xl-fast':         (1024, 20, _MK['xl']),
}

HELP = {
    'small-ultra':    "viruses, deepest: w4 c1 (all seeds)",
    'small-sensitive':"viruses, order+: w4 c2",   'small-normal':"viruses, family/genus: w4 c8",
    'small-fast':     "viruses, genus screen: w6 c12",
    'med-sensitive':  "bacteria/archaea: w5 c2",   'med-normal':"bacteria/archaea default: w6 c8",
    'med-fast':       "bacteria/archaea screen: w6 c30",
    'large-sensitive':"fungi/protists: w6 c2",     'large-normal':"fungi/protists: w6 c8",
    'large-fast':     "fungi/protists screen: w6 c20",
    'xl-sensitive':   "plant/animal (>0.5 Gaa): w6 c4",  'xl-normal':"plant/animal: w8 c8",
    'xl-fast':        "plant/animal screen: w7 c20",
}

DEFAULT_PRESET = 'med-normal'

MAX_PAIR_MATCHES = 32_000_000


def apply_preset(args):
    cfg = PRESETS.get(getattr(args, 'preset', None))
    if cfg is None:
        return

    args.w = cfg['w']
    args.shape = cfg['shape']

    if getattr(args, 'scaled', 10) == 10:
        args.scaled = cfg['scaled']

    _defaults = {
        'min_aai': 0.01,
        'max_evalue': 0.1,
        'min_chain_seeds': 2,
        'min_orf_cov': 0.20,
        'xdrop': 15,
        'max_pair_matches': 32_000_000,
    }
    for key, argparse_default in _defaults.items():
        current = getattr(args, key, None)
        if current is None or current == argparse_default:
            setattr(args, key, cfg[key])


def resolve_seed_params(args) -> Tuple[str, int]:
    shape: str = getattr(args, 'shape', None) or PRESETS[DEFAULT_PRESET]['shape']
    scaled: int = getattr(args, 'scaled', 10)
    return shape, scaled


def list_dbs(folder: str, ext: str = '.chaai') -> List[str]:
    if not os.path.isdir(folder):
        logging.error(f"GenomeDb folder not found: {folder}")
        sys.exit(1)
    return sorted(
        os.path.join(folder, f)
        for f in os.listdir(folder)
        if f.endswith(ext)
    )


def _mpm(args) -> int:
    v = getattr(args, 'max_pair_matches', None)
    return MAX_PAIR_MATCHES if v is None else v


def dist_options(args) -> dict:
    seed_shape, seed_scaled = resolve_seed_params(args)
    return dict(
        seed_shape=seed_shape, seed_scaled=seed_scaled,
        max_seed_occ=0, chain_band=25, chain_max_qgap=500,
        min_chain_seeds=getattr(args, 'min_chain_seeds', 2),
        min_aai=getattr(args, 'min_aai', 0.01),
        min_aligned_aa=0, min_chain_aligned_aa=0,
        min_matched_orfs=getattr(args, 'min_matched_orfs', 0),
        min_orf_cov=getattr(args, 'min_orf_cov', 0.20),
        aggressive_filter=getattr(args, 'aggressive_filter', False),
        xdrop=getattr(args, 'xdrop', 15),
        max_evalue=getattr(args, 'max_evalue', 0.1),
        max_pair_matches=_mpm(args),
        gene_scaled=getattr(args, 'gene_scaled', 1),
        temp_dir=getattr(args, 'temp_dir', '') or '',
    )
