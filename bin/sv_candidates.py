##########################
# BEGINING OF THE SCRIPT #
##########################

import sage.all
import sys
import os

if len(sys.argv) != 3:
    raise RuntimeError("usage: sv_candidates n d")
n = int(sys.argv[1])
d = int(sys.argv[2])

assert n > 0 and d > 0

def word_rank(u):
    return sum(2**i for i,j in enumerate(u) if j)

def word_unrank(n, length):
    u = []
    m = n
    while m:
        u.append(m%2)
        m //= 2
    assert len(u) <= length, "n = {} length = {} u ={}".format(n,length,u)
    u += [0]*(length-len(u))
    return tuple(int(i) for i in reversed(u))

def get_task(i, start, end, k):
    r"""
    Get the ``i``-th task between ``start`` and ``end``.
    """
    q = (end-start+1) // k
    r = (end-start+1) % k

    if i < r:
        return (start + (q+1)*i, start + (q+1)*(i+1) - 1)
    else:
        return (start + (q+1)*r + q*(i-r), start + (q+1)*r + q*(i-r+1) - 1)

nb_tasks = os.getenv('SLURM_NTASKS')
if nb_tasks is None:
    raise RuntimeError("SLURM_NTASKS not available")
nb_tasks = int(nb_tasks)

i = os.getenv('SLURM_PROCID')
if i is None:
    raise RuntimeError("SLURM_PROCID not available")
i = int(i)

if i < 0 or i >= nb_tasks:
    raise RuntimeError("i = {} should be in between 0 and nb_task={}".format(i, nb_tasks))

i_start, i_stop = get_task(i, 0, 2**n-1, nb_tasks)
u_start = word_unrank(i_start, n)
u_stop   = word_unrank(i_stop ,n)

from max_plus.sv_identities import sv_candidates

JOB_ID = os.getenv('SLURM_JOB_ID')
if JOB_ID is None:
    raise RuntimeError("no job id")

descr = '# sv candidates n={} d={} u_start={} u_stop={} #'.format(
        n, d, u_start, u_stop)
slurm  = '# JOB_ID: {}'.format(JOB_ID)
slurm += ' '*(len(descr)-len(slurm)-1) + '#\n'
t = '# NTASKS: {}'.format(nb_tasks)
slurm += t + ' '*(len(descr)-len(t)-1) + '#\n'
t = '# PROCID: {}'.format(i)
slurm += t + ' '*(len(descr)-len(t)-1) + '#'
print '#'*len(descr)
print descr
print slurm
print '#'*len(descr)
for u in sv_candidates(n, d, u_start, u_stop):
    print u
