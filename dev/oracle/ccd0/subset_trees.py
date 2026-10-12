#!/usr/bin/env python3
"""Subset a BEAST NEXUS tree file textually (newick lines kept byte-identical).
usage: subset_trees.py IN OUT K [burnin_fraction=0.1]
Keeps the header up to the first tree line, drops the first burnin fraction of
trees, then keeps K trees evenly spaced over the remainder (the last K if K
exceeds it), and closes the block with End;"""
import sys
inp, out, K = sys.argv[1], sys.argv[2], int(sys.argv[3])
burn = float(sys.argv[4]) if len(sys.argv) > 4 else 0.1
head, trees = [], []
with open(inp) as f:
    for line in f:
        s = line.lstrip().lower()
        if s.startswith("tree "):
            trees.append(line)
        elif not trees:
            head.append(line)
post = trees[int(len(trees) * burn):]
if K >= len(post):
    keep = post
else:
    step = len(post) / K
    keep = [post[int(i * step)] for i in range(K)]
with open(out, "w") as f:
    f.writelines(head)
    f.writelines(keep)
    f.write("End;\n")
print(inp, "trees", len(trees), "post-burnin", len(post), "kept", len(keep))
