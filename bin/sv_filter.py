import sys
import os
from max_plus import is_sv_identity

if len(sys.argv) != 3:
    raise ValueError
n = int(sys.argv[1])
d = int(sys.argv[2])

candidates = []
for filename in os.listdir('sv_relations/{}/{}'.format(d,n)):
    print "reading {}".format(filename)
    if os.path.isfile(filename):
        with open('filename') as f:
            lines = [line for line in f.read().split('\n') if line and not line.startsiwth('#')]
        candidates.append(map(eval,lines))

print '# n = {}  d = {}'.format(n,d)
print '# {} candidates'.format(len(candidates))
for t1,t2 in candidates:
    if is_sv_identity(t1,t2,d):
        print t1,t2

