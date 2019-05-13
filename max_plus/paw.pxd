#distutils: language=c++

from libcpp cimport bool as cppbool
from libcpp.string cimport string
from libcpp.vector cimport vector as cppvector
from libc.stdint cimport uint64_t

cdef extern from "<sstream>" namespace "std":
    cdef cppclass ostringstream:
        string str()

cdef extern from "paw_utils.h":
    ctypedef cppvector[uint64_t] int_vec
    uint64_t int_vec_sum(int_vec& u)
    void word_count(size_t& n0, size_t& n1, uint64_t w)
    size_t word_length(uint64_t)

    int cmp_lexico(int_vec&, int_vec&)
    int cmp_prod(int_vec&, int_vec&)

cdef extern from "paw.h":
    cppclass PAW:
        uint64_t w
        int_vec * n0
        int_vec * n1

        PAW()  except +
        PAW(uint64_t) except +
        PAW(uint64_t, int_vec&, int_vec&) except +

        size_t size()
        size_t length()

        void write()
        void set(uint64_t, int_vec&, int_vec&) except +

        cppbool operator==(const PAW&)
        cppbool operator!=(const PAW&)

    ostringstream& operator<<(ostringstream&, PAW&)


cdef class PartiallyAbelianizedWord(object):
    cdef PAW * paw


