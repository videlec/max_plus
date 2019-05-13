from sage.ext.memory_allocator cimport MemoryAllocator

from cpython.tuple cimport PyTuple_SET_ITEM, PyTuple_GET_ITEM
from cpython.int cimport PyInt_Check, PyInt_AS_LONG
from cpython.object cimport PyObject

def product_start_stop(tuple start, tuple stop):
    r"""
    Iterator over the words over {0,1} between start and stop inclusive.

    INPUT:

    - ``start``, ``stop`` -- the first and last lists

    EXAMPLES::

        sage: from max_plus.word import product_start_stop
        sage: for w in product_start_stop((0,0,0), (1,1,1)): print(w)
        (0, 0, 0)
        (0, 0, 1)
        (0, 1, 0)
        (0, 1, 1)
        (1, 0, 0)
        (1, 0, 1)
        (1, 1, 0)
        (1, 1, 1)
        sage: for w in product_start_stop((1,0,0,0), (1,1,0,1)): print(w)
        (1, 0, 0, 0)
        (1, 0, 0, 1)
        (1, 0, 1, 0)
        (1, 0, 1, 1)
        (1, 1, 0, 0)
        (1, 1, 0, 1)
    """
    cdef Py_ssize_t n = len(start)
    if len(stop) != n:
        raise ValueError("start, stop must have same length")
    
    if not all(PyInt_Check(i) for i in start):
        start = tuple(int(i) for i in start)
    if not all(PyInt_Check(i) for i in stop):
        stop = tuple(int(i) for i in stop)
    if not all(0 <= i <= 1 for i in start):
        raise ValueError
    if not all(0 <= i <= 1 for i in stop):
        raise ValueError

    if start > stop:
        return

    cdef tuple w = start[:]

    cdef Py_ssize_t j
    while w <= stop:
        yield w[:]
        
        j = n-1
        while j >= 0 and PyInt_AS_LONG(<object>PyTuple_GET_ITEM(w, j)) == 1:
            j -= 1
        if j == -1:
            return
        PyTuple_SET_ITEM(w, j, 1)
        j += 1
        while j < n:
            PyTuple_SET_ITEM(w, j, 0)
            j += 1
