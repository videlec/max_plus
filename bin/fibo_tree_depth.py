r"""
We compute the possible depths at which a given minimal
PAF appear word.
"""

from __future__ import print_function
from sage.all import *
from max_plus.paw import min_prod
from collections import defaultdict

for max_size in range(10,20):
    print("max_size = %d" % max_size)

    depth, M = min_prod([0,1], [0], max_size)

    f = open("depths%02d.txt" % max_size, "w")
    f.write("size:word:paf:depth\n")
    for i,m in enumerate(M[1:]):
        d = defaultdict(list)
        for p in m:
            d[p.word()].append(p)
        for w in sorted(d):
            for m in sorted(d[w], key=lambda x: depth[x]):
                f.write("%d:%d:%s:%s\n" % (i, depth[m], w, m))
    f.close()
