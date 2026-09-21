AA_LETTERS = "ARNDCQEGHILKMFPSTWYV*"
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = {"U": "C", "O": "K", "B": "N", "Z": "Q", "J": "L"}


def sanitize_protein(prot):
    if prot.endswith("*"):
        prot = prot[:-1]
    return "".join(AMBIGUOUS_AA.get(c, c) if (c in STANDARD_AA or c in AMBIGUOUS_AA) else "A"
                   for c in prot).replace("*", "A")


def unpack_5bit(packed, start, aa_len):
    out = []
    buf = bits = 0
    p = start
    for _ in range(aa_len):
        while bits < 5:
            buf |= packed[p] << bits; p += 1; bits += 8
        out.append(AA_LETTERS[buf & 0x1F]); buf >>= 5; bits -= 5
    return "".join(out)
