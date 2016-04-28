def parallel_unfold(args):
    r"""
    Helper for parallel functions.

    ``args`` should be a tuple ``(verbose, function, arguments)`` where

    - ``verbose`` -- a boolean indicating whether some additional

    - ``function`` -- the function to be parallelized
    """
    import sys
    import multiprocessing as mp
    from datetime import datetime
    from time import time

    verbose = args[0]
    f = args[1]
    args = args[2:]
    if verbose:
        t = datetime.now()
        t0 = time()
        print "{}: new job at {}:{}:{}\n  {}".format(
                mp.current_process().name, t.hour, t.minute,
                t.second, args)
        sys.stdout.flush()
    ans = f(*args)
    if verbose:
        print "{}: job done in {} seconds".format(mp.current_process().name, time()-t0)
        sys.stdout.flush()
    return ans

