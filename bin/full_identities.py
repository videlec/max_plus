import sys

import sage.all
from max_plus import *

d = int(sys.argv[1])
n = int(sys.argv[2])
print "d = {}  n = {}".format(d,n)

for u,v in full_identities_iterator(d,n,250,(1000 if d==2 else 50000),d==2):
    print "{} {}".format(''.join(map(str,u)), ''.join(map(str,v)))
