r"""
We compute the possible depths at which a given minimal
PAF appear word.
"""

from __future__ import print_function

import sage.all

import sys
import os

from max_plus.paw import min_prod
from collections import defaultdict

if len(sys.argv) == 2:
    nmin = int(sys.argv[1])
    nmax = nmin + 1
elif len(sys.argv) == 3:
    nmin = int(sys.argv[1])
    nmax = int(sys.argv[2])
else:
    raise ValueError("Usage: fibo_min_prod n1 [n2]")

for max_size in range(nmin, nmax):
    print("max_size = %d" % max_size)

    depth, fdepth, M = min_prod([0,1], [0], max_size)

    f = open("depths%02d.txt" % max_size, "w")
    f.write("size:length:depth:factor depth:word:paf\n")
    for i,m in enumerate(M[1:]):
        d = defaultdict(list)
        for p in m:
            d[p.word()].append(p)
        for w in sorted(d):
            for m in sorted(d[w], key=lambda x: depth[x]):
                f.write("%d:%d:%d:%d:%s:%s\n" % (i, m.length(), depth[m], fdepth[m], w, m))
    f.close()
