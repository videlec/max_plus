import time

def write_line(f, s):
    if s is None:
        f.write('#'*73)
        f.write('\n')
    else:
        f.write('# ')
        f.write(s)
        f.write(' '*(70-len(s)))
        f.write('#\n')

def get_task(i, start, end, k):
    r"""
    Get the ``i``-th task between ``start`` and ``end``.
    """
    i = int(i)
    start = int(start)
    end = int(end)
    k = int(k)

    q = (end-start+1) // k
    r = (end-start+1) % k

    if end < start:
        raise ValueError("end (={}) is smaller or equal than start (={})".format(end,start))

    if i < r:
        out = (start + (q+1)*i, start + (q+1)*(i+1) - 1)
    else:
        out = (start + (q+1)*r + q*(i-r), start + (q+1)*r + q*(i-r+1) - 1)

    return out

def word_nice_str(u):
    return ''.join(map(str,u))

def word_rank(u):
    return sum(2**i for i,j in enumerate(u) if j)

def word_unrank(n, length):
    n = int(n)
    length = int(length)
    u = []
    m = n
    while m:
        u.append(m%2)
        m //= 2
    assert len(u) <= length, "n = {} length = {} u ={}".format(n,length,u)
    u += [0]*(length-len(u))
    return tuple(int(i) for i in reversed(u))

def run_task(arg):
    r"""
    Run *one* task.

    This is the only part that actually loads Sage.
    """
    typ, n, d, (i_start, i_stop) = arg

    import multiprocessing as mp
    import sage.all
    if typ == 'sv':
        from max_plus.sv_identities import sv_candidates as candidates
    elif typ == 'vv':
        from max_plus.vv_identities import vv_candidates as candidates
    else:
        raise ValueError("typ should either be sv or vv")

    return list(candidates(n, d, word_unrank(i_start,n), word_unrank(i_stop,n)))

MPITAG_MASTER_SUBMIT_JOB=12
MPITAG_SLAVE_RETURN_RESULT=14
MPITAG_SLAVE_AVAILABLE=13

MPI_CLOSE_SLAVE=None

class Worker:
    def __init__(self, comm, typ, n, d):
        import multiprocessing as mp
        self.comm = comm
        self.ncpus = mp.cpu_count()
        self.typ = typ
        self.n = n
        self.d = d

    def main(self):
        data = self.comm.recv(source=0, tag=MPITAG_MASTER_SUBMIT_JOB)
        while data != MPI_CLOSE_SLAVE:
            start,stop = data
            self.run_parallel(start, stop)
            self.comm.send(None, dest=0, tag=MPITAG_SLAVE_AVAILABLE)
            data = self.comm.recv(source=0, tag=MPITAG_MASTER_SUBMIT_JOB)

    def run_parallel(self, i_start, i_stop):
        import multiprocessing as mp

        u_start = word_unrank(i_start, self.n)
        u_stop = word_unrank(i_stop, self.n)

        s_start = word_nice_str(u_start)
        s_stop = word_nice_str(u_stop)

        outdir = os.path.join('{}_candidates'.format(self.typ), str(self.d), str(self.n))

        filename = os.path.join(outdir, '{}-{}'.format(s_start, s_stop))
        if os.path.isfile(filename):
            f = open(filename)
            done = 'END' in f.read()
            f.close()
            if done:
                return

        f = open(filename, 'w')
        write_line(f, None)
        write_line(f, '{} candidates n={} d={}'.format(self.typ,self.n,self.d))
        write_line(f, 'u_start: {} ({})'.format(s_start, i_start))
        write_line(f, 'u_stop : {} ({})'.format(s_stop, i_stop))
        write_line(f, None)

        nb_subtasks = (i_stop-i_start+1) // 5
        tasks = [(self.typ, self.n, self.d, get_task(j, i_start, i_stop, nb_subtasks)) for j in xrange(nb_subtasks)]
        pool = mp.Pool(self.ncpus)
        for res in pool.imap_unordered(run_task, tasks):
            for u0,u1 in res:
                f.write('{} {}\n'.format(u0, u1))
            f.flush()
        write_line(f, None)
        write_line(f, 'END')
        write_line(f, None)
        f.close()

class Scheduler:
    def __init__(self, comm, job_list):
        self.comm = comm
        from mpi4py import MPI
        self.rank = self.comm.Get_rank()
        self.size = self.comm.Get_size()

        self.submitted = [None]*self.size
        self.job_list  = job_list
        self.results   = {}

    def submit_job(self, rank=None):
        if rank is None:
            for k in range(1,self.size):
                self.submit_job(k)

        elif self.job_list:
            job = self.job_list.pop()
            if self.submitted[rank] is not None:
                raise RuntimeError
            self.submitted[rank] = (job, time.time())
            self.comm.isend(job, dest=rank, tag=MPITAG_MASTER_SUBMIT_JOB)

    def gather_results(self, verbose=False):
        for k in range(1,self.size):
            job = self.submitted[k]
            if job is None:
                continue
            job,t0 = job
            if self.comm.Iprobe(source=k, tag=MPITAG_SLAVE_AVAILABLE):
                print "worker {} finished task {} in {} secs".format(k, job, time.time()-t0)
                self.comm.recv(source=k, tag=MPITAG_SLAVE_AVAILABLE)
                self.submitted[k] = None
                if self.job_list:
                    self.submit_job(k)
                else:
                    self.close(k)

    def close(self, rank):
        self.comm.isend(MPI_CLOSE_SLAVE, dest=rank,
                tag=MPITAG_MASTER_SUBMIT_JOB)

    def done(self):
        return not self.job_list and all(x is None for x in self.submitted)

    def get_results(self):
        return self.results

if __name__ == '__main__':
    import sys
    import os
    from mpi4py import MPI

    if len(sys.argv) != 4:
        raise RuntimeError("usage: candidates sv|vv n d")
    typ = sys.argv[1]
    n = int(sys.argv[2])
    d = int(sys.argv[3])

    # tags
    comm = MPI.COMM_WORLD

    if comm.Get_rank() == 0:
        # the master process
        i_start = 2**(n-1)
        i_stop = 2**n-1
        nb_tasks = min(comm.Get_size()**3, 2**(n-1))
        joblist = [get_task(k, i_start, i_stop, nb_tasks) for k in range(nb_tasks)]
        scheduler = Scheduler(comm, joblist)
        scheduler.submit_job()

        while not scheduler.done():
            scheduler.gather_results(verbose=True)
            time.sleep(0.5)

    else:
        # slave processes
        worker = Worker(comm, typ, n, d)
        worker.main()
