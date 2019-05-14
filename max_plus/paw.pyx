# distutils: language=c++
# distutils: sources=max_plus/paw_utils.cpp max_plus/paw_cpp.cpp
r"""
Compute the minimal products
"""

from cpython.object cimport Py_LT, Py_LE, Py_EQ, Py_NE, Py_GT, Py_GE

from libc.stdint cimport uint64_t

from libcpp.list cimport list as cpplist

from itertools import combinations, product
from sage.combinat.integer_lists.invlex import IntegerListsLex
from sage.combinat.words.all import FiniteWords

W = FiniteWords([0,1])
WordChar = W._element_classes['char']

cdef int int_vec_from_python(int_vec* v, l) except -1:
    cdef size_t n = len(l)
    cdef size_t i
    cdef uint64_t j
    v.resize(n)
    for i,j in enumerate(l):
        v[0][i] = j
    return 0

cdef class PartiallyAbelianizedWord(object):
    r"""
    A partially abelianzed word.

    EXAMPLES::

        sage: from max_plus.paw import PartiallyAbelianizedWord

        sage: PartiallyAbelianizedWord([0,1], [0], [0])
        01
        sage: PartiallyAbelianizedWord([0,1], [1], [0])
        0(1,0)1
        sage: PartiallyAbelianizedWord([0,1], [0], [1])
        0(0,1)1
        sage: PartiallyAbelianizedWord([0,1], [4], [1])
        0(4,1)1

        sage: PartiallyAbelianizedWord([0,1,1,0], [1,2,3], [3,3,1])
        0(1,3)1(2,3)1(3,1)0
    """
    def __cinit__(self):
        self.paw = new PAW();

    def __dealloc__(self):
        del self.paw

    def __hash__(self):
        return self.paw.hash()

    def __init__(self, w, ab0=None, ab1=None):
        r"""
        TESTS::

            sage: from max_plus.paw import PartiallyAbelianizedWord
            sage: PartiallyAbelianizedWord("00")
            00
            sage: PartiallyAbelianizedWord("0(0,0)0")
            00
            sage: PartiallyAbelianizedWord("0(1,1)001(1,2)0")
            0(1,1)001(1,2)0
        """
        if ab0 is None or ab1 is None:
            if not isinstance(w, str):
                raise TypeError("must be a string")
            s = w
            w = []
            ab0 = []
            ab1 = []
            m = 0
            while m < len(s):
                k1 = s.find('(', m)
                if k1 == -1:
                    w.extend(map(int, s[m:]))
                    ab0.extend([0] * (len(s)-m-1))
                    ab1.extend([0] * (len(s)-m-1))
                    break
                else:
                    k2 = s.find(',', k1)
                    k3 = s.find(')', k2)
                    w.extend(map(int, s[m:k1]))
                    ab0.extend([0] * (k1-m-1))
                    ab1.extend([0] * (k1-m-1))
                    ab0.append(int(s[k1+1:k2]))
                    ab1.append(int(s[k2+1:k3]))
                    m = k3 + 1

        assert len(w) -1 == len(ab0) == len(ab1)

        cdef uint64_t u,i,j
        u = 1
        for j in reversed(w):
            if j != 0 and j != 1:
                raise ValueError
            u <<= 1
            u += j

        cdef int_vec * n0 = new int_vec()
        cdef int_vec * n1 = new int_vec()

        int_vec_from_python(n0, ab0)
        int_vec_from_python(n1, ab1)

        self.paw.set(u, n0[0], n1[0])

        del n0
        del n1

    def reverse(self):
        r"""
        EXAMPLES::

            sage: from max_plus.paw import PartiallyAbelianizedWord
            sage: w = PartiallyAbelianizedWord([1,0,0,1,0], [0,1,2,3],[3,3,1,1])
            sage: w.reverse()
            0(3,1)1(2,1)0(1,3)0(0,3)1
        """
        cdef size_t i
        r = self.word()
        n0 = [self.paw.n0[0][i] for i in range(self.paw.n0.size())]
        n1 = [self.paw.n1[0][i] for i in range(self.paw.n1.size())]
        return PartiallyAbelianizedWord(r[::-1], n0[::-1], n1[::-1])

    def word(self):
        r"""
        EXAMPLES::

            sage: from max_plus.paw import PartiallyAbelianizedWord
            sage: PartiallyAbelianizedWord([0,1,1,0], [1,2,3], [3,3,1]).word()
            word: 0110
        """
        cdef list l = []
        cdef uint64_t u = self.paw.w
        while u != 1:
            l.append(u&1)
            u >>= 1
        return WordChar(W, l)

    def gaps(self):
        r"""
        EXAMPLES::

            sage: from max_plus.paw import PartiallyAbelianizedWord
            sage: PartiallyAbelianizedWord([0,1,1,0], [1,2,3], [3,3,1]).gaps()
            [(1, 3), (2, 3), (3, 1)]
        """
        from sage.modules.all import FreeModule
        from sage.rings.integer_ring import ZZ
        V = FreeModule(ZZ, 2)

        ret = []
        cdef size_t i
        for i in range(self.paw.n0.size()):
            ret.append(V((self.paw.n0[0][i], self.paw.n1[0][i])))
        return ret

    def size(self):
        r"""
        Return the number of marked positions.

        EXAMPLES::

            sage: from max_plus.paw import PartiallyAbelianizedWord
            sage: PartiallyAbelianizedWord([0,1,1,0], [1,2,3], [3,3,1]).size()
            4
        """
        return self.paw.size()

    def length(self):
        r"""
        Return the total length.

        EXAMPLES::

            sage: from max_plus.paw import PartiallyAbelianizedWord
            sage: PartiallyAbelianizedWord([0,1,1,0], [1,2,3], [3,3,1]).length()
            17
        """
        return self.paw.length()

    def __richcmp__(self, other, op):
        r"""
        TESTS::

            sage: from max_plus.paw import PartiallyAbelianizedWord as PAW
            sage: PAW([0], [], []) == PAW([0], [], [])
            True
            sage: PAW([0], [], []) == PAW([1], [], [])
            False

            sage: PAW([0,0], [0], [1]) == PAW([0,0], [0], [1])
            True
            sage: PAW([0,0], [0], [1]) == PAW([0,0], [1], [0])
            False
            sage: PAW([0,0], [0], [1]) == PAW([0,0], [0], [0])
            False
            sage: PAW([0,0], [1], [0]) == PAW([0,1], [1], [0])
            False

            sage: PAW([0,0], [0], [1]) != PAW([0,0], [0], [1])
            False
            sage: PAW([0,0], [0], [1]) != PAW([0,0], [1], [0])
            True
            sage: PAW([0,0], [0], [1]) != PAW([0,0], [0], [0])
            True

        """
        if type(self) is not type(other) or (op != Py_EQ and op != Py_NE):
            raise TypeError

        if op == Py_EQ:
            return (<PartiallyAbelianizedWord> self).paw[0] == (<PartiallyAbelianizedWord> other).paw[0]
        else:
            return (<PartiallyAbelianizedWord> self).paw[0] != (<PartiallyAbelianizedWord> other).paw[0]

    def __repr__(self):
        cdef ostringstream s
        s << self.paw[0]
        return s.str()

    def abelian_vector(self):
        r"""
        Return the Abelian vector.

        EXAMPLES::

            sage: from max_plus.paw import PartiallyAbelianizedWord as PAW
            sage: PAW([0], [], []).abelian_vector()
            (1L, 0L)
            sage: PAW([0,1], [3], [7]).abelian_vector()
            (4L, 8L)
            sage: PAW([0,0,0], [2,2], [0,1]).abelian_vector()
            (7L, 1L)
        """
        cdef size_t n0 = 0
        cdef size_t n1 = 0
        word_count(n0, n1, self.paw.w)
        return (n0 + int_vec_sum(self.paw.n0[0]), n1 + int_vec_sum(self.paw.n1[0]))

    def match(self, list w):
        r"""
        Check whether this partial abelianization comes from ``w``

        EXAMPLES::

            sage: from max_plus.paw import PartiallyAbelianizedWord
            sage: from itertools import product
            sage: p = PartiallyAbelianizedWord([1,1],[0],[0])
            sage: for l in range(10):
            ....:     for u in product([0,1], repeat=l):
            ....:         if p.match(list(u)):
            ....:             print(u)
            (1, 1)

            sage: p = PartiallyAbelianizedWord([0,0],[1],[1])
            sage: for l in range(10):
            ....:     for u in product([0,1], repeat=l):
            ....:         if p.match(list(u)):
            ....:             print(u)
            (0, 0, 1, 0)
            (0, 1, 0, 0)

            sage: p = PartiallyAbelianizedWord([0,0,0],[1,1],[1,2])
            sage: for l in range(10):
            ....:     for u in product([0,1], repeat=l):
            ....:         if p.match(list(u)):
            ....:             print(u)
            (0, 0, 1, 0, 0, 1, 1, 0)
            (0, 0, 1, 0, 1, 0, 1, 0)
            (0, 0, 1, 0, 1, 1, 0, 0)
            (0, 1, 0, 0, 0, 1, 1, 0)
            (0, 1, 0, 0, 1, 0, 1, 0)
            (0, 1, 0, 0, 1, 1, 0, 0)
        """
        if len(w) != self.paw.length() or self.paw.w&1 != w[0]:
            return False

        cdef size_t i
        cdef int l
        cdef uint64_t u = self.paw.w
        cdef list s
        u >>= 1
        cp = 1 # current position
        for i in range(self.paw.n0.size()):
            l = self.paw.n0[0][i] + self.paw.n1[0][i]
            s = w[cp:cp+l]
            if (s.count(0) != self.paw.n0[0][i] or
                s.count(1) != self.paw.n1[0][i] or
                w[cp+l] != u&1):
                return False
            cp += l+1
            u >>= 1
        assert u == 1
        return True

    def substitutions(self, w0, w1, max_size=None):
        r"""
        Iterate through the substitutions of this PAW under the morphism
        0 -> w0 and 1 -> w1.

        EXAMPLES::

            sage: from max_plus.paw import PartiallyAbelianizedWord

            sage: w = PartiallyAbelianizedWord([0], [], [])
            sage: for u in w.substitutions([0,1],[0]):
            ....:     print(u)
            ....:
            01
            0
            1

            sage: w = PartiallyAbelianizedWord([1], [], [])
            sage: for u in w.substitutions([0,1],[0]):
            ....:     print(u)
            0

            sage: w = PartiallyAbelianizedWord([0,0], [0], [0])
            sage: for u in w.substitutions([0,1],[0]):
            ....:     print(u)
            0101
            010
            01(1,0)1
            0(0,1)01
            101
            0(0,1)0
            0(1,1)1
            10
            1(1,0)1

            sage: w = PartiallyAbelianizedWord([0,0], [1], [1])
            sage: list(w.substitutions([0,1],[0], max_size=2)) == [u for u in w.substitutions([0,1],[0]) if u.size() <= 2]
            True

            sage: w = PartiallyAbelianizedWord([0,1,0], [0,0], [0,0])
            sage: for u in w.substitutions([0,1],[0]):
            ....:     print(u)
            01001
            0100
            010(1,0)1
            0(0,1)001
            1001
            0(0,1)00
            0(0,1)0(1,0)1
            100
            10(1,0)1

            sage: w = PartiallyAbelianizedWord([0,0], [1],[1])
            sage: L1 = list(w.substitutions([0,1],[0], max_size=2))
            sage: L2 = [u for u in w.substitutions([0,1],[0]) if u.size() <= 2]
            sage: L1 == L2
            True
        """
        w0 = WordChar(W, w0)
        w1 = WordChar(W, w1)

        if not w0 or not w1:
            raise ValueError("images can not be empty")
        cdef Py_ssize_t l0 = len(w0)
        cdef Py_ssize_t l1 = len(w1)

        # ARGH: it is much better to take a word from sage.combinat.words
        # since we pick a lot of slices below!!
        sw = self.word()
        if max_size is None:
            max_size = sum(l1 if letter else l0 for letter in sw)

        # substitution matrix
        cdef int a = w0.count(0)
        cdef int c = w0.count(1)
        cdef int b = w1.count(0)
        cdef int d = w1.count(1)

        cdef list sub = [w1 if letter else w0 for letter in sw]
#        print('sub = {}'.format(sub))

        cdef int j, n0_left, n1_left
        cdef list w
        cdef tuple p
        cdef PartiallyAbelianizedWord r

        # now we choose all possible combination of letters among the images
        for L in IntegerListsLex(
            length=len(sw),
            min_part=1,
            max_sum=max_size,
            ceiling=[len(u) for u in sub]):

            for p in product(*[combinations(range(len(s)),l) for s,l in zip(sub, L)]):
#                print('new sub with p = {}'.format(p))
                # treat appart the first
                r = PartiallyAbelianizedWord.__new__(PartiallyAbelianizedWord)

                pi = p[0]
                subi = sub[0]
                w = [subi[j] for j in pi]

                for j in range(len(pi)-1):
                    r.paw.n0.push_back(subi[pi[j]+1:pi[j+1]].count(0))
                    r.paw.n1.push_back(subi[pi[j]+1:pi[j+1]].count(1))

                n0_left = subi[pi[-1]+1:].count(0)
                n1_left = subi[pi[-1]+1:].count(1)

                for i in range(self.paw.n0.size()):
                    pi = p[i+1]
                    subi = sub[i+1]
#                    print(' i = {}'.format(i))
#                    print(' pi = {}'.format(pi))
#                    print(' sub[i] = {}'.format(sub[i+1]))
#                    print(' n_left = {}'.format(n_left))
                    w.extend(subi[j] for j in pi)
                    # a * n0[i] + b * n1[i]
                    suff = subi[:pi[0]]
                    # n.append(n_left + mat * gaps[i] + V(subi[:pi[0]].abelian_vector()))

                    r.paw.n0.push_back(n0_left + a * self.paw.n0[0][i] + b * self.paw.n1[0][i] + suff.count(0))
                    r.paw.n1.push_back(n1_left + c * self.paw.n0[0][i] + d * self.paw.n1[0][i] + suff.count(1))
                    for j in range(len(pi)-1):
                        r.paw.n0.push_back(subi[pi[j]+1:pi[j+1]].count(0))
                        r.paw.n1.push_back(subi[pi[j]+1:pi[j+1]].count(1))

                    n0_left = subi[pi[-1]+1:].count(0)
                    n1_left = subi[pi[-1]+1:].count(1)

                r.paw.w = 1
                for j in reversed(w):
                    r.paw.w <<= 1
                    r.paw.w += j

                yield r

cpdef clean_update(list L, PartiallyAbelianizedWord v):
    r"""
    Remove all elements in L that are larger than v (for the product order)
    and add v to L if it does not have anybody smaller.

    OUTPUT: new list or not

    EXAMPLES::

        sage: from max_plus.paw import PartiallyAbelianizedWord, clean_update

        sage: L = [PartiallyAbelianizedWord([0,0], [1], [1])]
        sage: clean_update(L, PartiallyAbelianizedWord([0,0], [1], [1]))
        sage: clean_update(L, PartiallyAbelianizedWord([0,0], [2], [1]))
        sage: L
        [0(1,1)0]
        sage: clean_update(L, PartiallyAbelianizedWord([0,0], [1], [0]))
        sage: clean_update(L, PartiallyAbelianizedWord([0,0], [0], [1]))
        sage: clean_update(L, PartiallyAbelianizedWord([0,0], [1], [1]))
        sage: L
        [0(1,0)0, 0(0,1)0]
        sage: clean_update(L, PartiallyAbelianizedWord([0,0], [0], [0]))
        sage: L
        [00]

        sage: L = [PartiallyAbelianizedWord([0,0,0],[0,0],[1,0])]
        sage: clean_update(L, PartiallyAbelianizedWord([0,0,0], [0,0], [0,1]))
        sage: L
        [0(0,1)00, 00(0,1)0]
    """
    cdef Py_ssize_t i = 0
    cdef PartiallyAbelianizedWord u
    cdef int test0, test1
    while i < len(L):
        u = L[i]
        assert u.paw.w == v.paw.w
        test0 = cmp_prod(u.paw.n0[0], v.paw.n0[0])
        test1 = cmp_prod(u.paw.n1[0], v.paw.n1[0])
        if test0 == 2 or test1 == 2:
            # incomparable
            i += 1
        elif test0 <= 0 and test1 <= 0:
            # found somebody smaller or equal
            return
        elif test0 >= 0 and test1 >= 0:
            # found somebody larger
            del L[i]
        else:
            # incomparable
            i += 1

    L.append(v)


cpdef bint lt_size(PartiallyAbelianizedWord w1, PartiallyAbelianizedWord w2):
    # compare sizes
    cdef size_t s1 = w1.paw.size()
    cdef size_t s2 = w2.paw.size()
    if s1 < s2:
        return True
    elif s1 > s2:
        return False

    # compare lengths
    cdef size_t l1 = w1.paw.length()
    cdef size_t l2 = w2.paw.length()
    if l1 < l2:
        return True
    if l1 > l2:
        return False

    # don't care
    return False

cdef uint64_t all_one(uint64_t u):
    cdef uint64_t v = 0
    while u:
        v <<= 1
        v += 1
        u >>= 1
    return v

cdef void argmin_size(list L, Py_ssize_t m, Py_ssize_t * m_min, Py_ssize_t * i_min):
    r"""
    A minimum element for comparison by lengths and sizes.
    """
    # 
    cdef Py_ssize_t lim = all_one(m)
    cdef Py_ssize_t i = 0
    cdef PartiallyAbelianizedWord u = L[m][0]
    cdef PartiallyAbelianizedWord v
    cdef list s
    m_min[0] = m
    i_min[0] = 0
    while m <= lim:
        s = L[m]
        for i in range(len(s)):
            v = s[i]
            if lt_size(v, u) == True:
                m_min[0] = m
                i_min[0] = i
                u = v
        m += 1

cpdef Py_ssize_t find_smaller_prod(list L, PartiallyAbelianizedWord v):
    r"""
    Find an element smaller than ``v`` (for the product order) in the list ``L``.

    Return the index of a smaller element if any or ``-1``.

    EXAMPLES::

        sage: from max_plus.paw import PartiallyAbelianizedWord as PAW, find_smaller_prod
        sage: u0 = PAW([0]*3, [1,1], [1,3])
        sage: u1 = PAW([0]*3, [3,1], [1,1])
        sage: u2 = PAW([0]*3, [1,2], [1,2])

        sage: v1 = PAW([0]*3, [1,1], [2,3])
        sage: v2 = PAW([0]*3, [3,1], [1,1])
        sage: v3 = PAW([0]*3, [2,2], [2,2])
        sage: v4 = PAW([0]*3, [1,1], [1,1])

        sage: L = [u0, u1, u2]
        sage: find_smaller_prod(L, v1)
        0
        sage: find_smaller_prod(L, v2)
        1
        sage: find_smaller_prod(L, v3)
        2
        sage: find_smaller_prod(L, v4)
        -1
    """
    cdef Py_ssize_t i
    cdef PartiallyAbelianizedWord u
    for i in range(len(L)):
        u = L[i]
        if cmp_prod(u.paw.n0[0], v.paw.n0[0]) <= 0 and cmp_prod(u.paw.n1[0], v.paw.n1[0]) <= 0:
            return i
    return -1

def min_prod(s0, s1, max_size):
    r"""
    Return the set of minimal elements for the product order in the language
    of the substitution ``0 -> s0`` and ``1 -> s1``.

    EXAMPLES:

    The extrema for the Fibonacci language::

        sage: from max_plus.paw import min_prod
        sage: M = min_prod([0,1], [0], 5)
        sage: M[1]
        [0, 1]
        sage: M[2]
        [00, 10, 01, 1(1,0)1]
        sage: M[3]
        [0(0,1)00,
         00(0,1)0,
         100,
         010,
         1(1,0)10,
         001,
         101,
         01(1,0)1,
         1(2,0)1(1,0)1,
         1(1,0)1(2,0)1]

    Checking invariance with respect to the mirror::

        sage: for k in range(2,5):
        ....:     for w in M[k]:
        ....:         assert w.reverse() in M[k]
    """
    # at each time of the loops all the elements in todo and mins are
    # incomparable
    # lt_size provides a total order: we would better use an ordered
    # queue for "todo"... but then it will not be that simple
    # to apply clean update

    cdef list mins = [[] for _ in range(2**(max_size+1))]
    cdef list todo = [[] for _ in range(2**(max_size+1))]
    cdef Py_ssize_t m = 2
    cdef PartiallyAbelianizedWord u, v
    cdef Py_ssize_t n = 0
    cdef Py_ssize_t i = 0

    cdef dict depth = {}
    root = PartiallyAbelianizedWord([0],[],[])
    todo[2].append(root)
    depth[root] = 0

    while True:
        while m < len(todo) and not todo[m]:
            m += 1
        if m == len(todo):
            break

        # pick the minimum in size
        argmin_size(todo, m, &n, &i)
        u = todo[n].pop(i)

#        print('mins: {}'.format(mins))
#        print('todo: {}'.format(todo))
#        print('u   : {}'.format(u))

        # safety check
#        assert find_smaller_prod(mins[n], u) == -1

        mins[n].append(u)

        # substitute
        for v in u.substitutions(s0, s1, max_size=max_size):
            if v not in depth:
                depth[v] = depth[u] + 1
            else:
                depth[v] = min(depth[v], depth[u]+1)
            if v.paw.w < m:
                m = v.paw.w
            if find_smaller_prod(mins[v.paw.w], v) == -1:
                clean_update(todo[v.paw.w], v)

        # safety check
#        assert u not in todo[u.paw.w]

    # indices 2,3 -> size 1
    # indices 4,5,6,7 -> size 2
    # etc
    M = [None]
    j = 2
    for i in range(1, max_size+1):
        M.append([])
        for k in range(2**i):
            M[-1].extend(mins[j+k])
        j+= 2**i

    return depth, M
