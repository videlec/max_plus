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
files = os.listdir(indir)
print "reading {} files...".format(len(files)),
n_errors = 0
not_finished = 0
for filename in files:
    filename = os.path.join(indir,filename)
    if os.path.isfile(filename):
        not_end = True
        with open(filename) as f:
            line = f.readline()
            while line:
                if not line.startswith('#'):
                    line = line.split()
                    if len(line) != 2 or len(line[0]) != n or len(line[1]) != n:
                        n_errors += 1
                    else:
                        t0,t1 = line
                        candidates.append((map(int,t0), map(int,t1)))
                elif 'END' in line:
                    not_end = False
                line = f.readline()
            not_finished += not_end
print "done"
candidates.sort()
print '{} candidates'.format(len(candidates))
if n_errors:
    print "-> {} malformated lines".format(n_errors)
if not_finished:
    print "-> {} computation not terminated".format(not_finished)
if n_errors or not_finished:
    output = None
else:
    output = open(os.path.join('sv_identities/{}_{}'.format(d,n)), 'w')

num=0
for t1,t2 in candidates:
    if is_sv_identity(t1,t2,d):
        rel = '{} = {}'.format(''.join(map(str,t1)),''.join(map(str,t2)))
        if output is not None:
            output.write(rel + '\n')
        print rel
        num += 1
print '{} identities'.format(num)
if output is not None:
    output.close()
