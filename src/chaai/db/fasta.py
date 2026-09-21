def read_fasta(path):
    seqs, name, cur = [], None, []
    for line in open(path):
        if line.startswith(">"):
            if name is not None:
                seqs.append((name, "".join(cur)))
            name = line[1:].split()[0]; cur = []
        else:
            cur.append(line.strip())
    if name is not None:
        seqs.append((name, "".join(cur)))
    return seqs
