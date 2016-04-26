import sage.all
import sys
import os
from max_plus import is_sv_identity

if len(sys.argv) != 3:
    raise ValueError
n = int(sys.argv[1])
d = int(sys.argv[2])

indir = 'sv_candidates/{}/{}'.format(d,n)

candidates = []
print "reading files..."
for filename in os.listdir(indir):
    filename = os.path.join(indir,filename)
    if os.path.isfile(filename):
        with open(filename) as f:
            line = f.readline()
            while line:
                if not line.startswith('#'):
                    t0,t1 = line.split()
                    assert len(t0) == len(t1) == n, (t0, t1)
                    candidates.append((map(int,t0), map(int,t1)))
                line = f.readline()
candidates.sort()
print "done"
print 'n = {}  d = {}'.format(n,d)
print '{} candidates'.format(len(candidates))
num=0
output = open(os.path.join('sv_identities/{}_{}'.format(n,d)), 'w')
for t1,t2 in candidates:
    if is_sv_identity(t1,t2,d):
        output.write('{} {}\n'.format(''.join(map(str,t1)),''.join(map(str,t2))))
        num += 1
print 'got {} identities'.format(num)
output.close()
