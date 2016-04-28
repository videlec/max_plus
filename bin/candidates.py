def word_nice_str(u):
    return ''.join(map(str,u))

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

    if end < start:
        raise ValueError("end (={}) is smaller or equal than start (={})".format(end,start))

    if i < r:
        out = (start + (q+1)*i, start + (q+1)*(i+1) - 1)
    else:
        out = (start + (q+1)*r + q*(i-r), start + (q+1)*r + q*(i-r+1) - 1)

    return out

def write_line(f, s):
    if s is None:
        f.write('#'*73)
        f.write('\n')
    else:
        f.write('# ')
        f.write(s)
        f.write(' '*(70-len(s)))
        f.write('#\n')

def run_task(arg):
    typ, n, d, (i_start, i_stop), outdir, jobid = arg

    import sage.all
    if typ == 'sv':
        from max_plus.sv_identities import sv_candidates as candidates
    elif typ == 'vv':
        from max_plus.vv_identities import vv_candidates as candidates
    else:
        raise ValueError("typ should either be sv or vv")

    u_start = word_unrank(i_start, n)
    u_stop = word_unrank(i_stop, n)

    p = mp.current_process()
    s_start = word_nice_str(u_start)
    s_stop = word_nice_str(u_stop)
    filename = os.path.join(outdir, '{}-{}'.format(s_start, s_stop))
    if os.path.isfile(filename):
        f = open(filename)
        done = 'END' in f.read()
        f.close()
        if done:
            return

    f = open(filename, 'w')
    write_line(f, None)
    write_line(f, '{} candidates n={} d={}'.format(typ,n,d))
    write_line(f, 'u_start: {} ({})'.format(s_start, i_start))
    write_line(f, 'u_stop : {} ({})'.format(s_stop, i_stop))
    write_line(f, 'JOB_ID : {}'.format(jobid))
    write_line(f, 'NTASKS : {}'.format(nb_tasks))
    write_line(f, 'PROCID : {}'.format(i))
    write_line(f, p.name)
    write_line(f, None)
    for u in candidates(n, d, u_start, u_stop):
        f.write('{} {}\n'.format(word_nice_str(u[0]), word_nice_str(u[1])))
        f.flush()
    write_line(f, None)
    write_line(f, 'END')
    write_line(f, None)
    f.close()

if __name__ == '__main__':
    import sys
    import os
    import multiprocessing as mp
    from time import time

    if len(sys.argv) != 4 and len(sys.argv) != 5:
        raise RuntimeError("usage: candidates sv|vv n d [MAX_IDENTITIES_PER_TASK]")
    typ = sys.argv[1]
    n = int(sys.argv[2])
    d = int(sys.argv[3])
    if len(sys.argv) == 5:
        max_id = int(sys.argv[4])
    else:
        max_id = float('inf')

    assert n > 0 and d > 0 and max_id >= 1

    nb_tasks = os.getenv('SLURM_NTASKS')
    if nb_tasks is None:
        raise RuntimeError("SLURM_NTASKS not available")
    nb_tasks = int(nb_tasks)

    i = os.getenv('SLURM_PROCID')
    if i is None:
        raise RuntimeError("SLURM_PROCID not available")
    i = int(i)

    jobid = os.getenv('SLURM_JOB_ID')
    if jobid is None:
        raise RuntimeError("SLURM_JOB_ID not available")

    if i < 0 or i >= nb_tasks:
        raise RuntimeError("i = {} should be in between 0 and nb_task={}".format(i, nb_tasks))

    ncpus = mp.cpu_count()

    outdir = os.path.join('{}_candidates'.format(typ), str(d), str(n))

    # (i_start, i_stop) is the interval allocated to this task
    # we further divide it according to the number of cpus available
    i_start, i_stop = get_task(i, 2**(n-1), 2**n-1, nb_tasks)

    tasks = []
    nb_subtasks = max(min(ncpus**2, i_stop-i_start+1), int((i_stop - i_start + 1)/max_id))
    print "TASK {} cut into {} subtasks (from {} to {})".format(i, nb_subtasks, i_start, i_stop)
    tasks = ((typ, n, d, get_task(j, i_start, i_stop, nb_subtasks), outdir, jobid) \
             for j in xrange(nb_subtasks))

    t0 = time()
    pool = mp.Pool(ncpus)
    for _ in pool.imap_unordered(run_task, tasks):
        pass
    pool.terminate()
    pool.join()
    t0 = time() - t0

    print "TASK {} done in {} secs".format(i, t0)
