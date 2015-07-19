"""
Some experiments on symbolic max plus matrices

We consider semigroup relations for various max-plus matrix semigroups:

 * upper triangular matrices with identical diagonals ('tri_sim_diag')
 * upper triangular matrices ('tri')
 * all matrices ('all')

If u = v is a relation then

pus = pvs is also a relation...

So we should look for primitive relations.


1. For band 2x2 matrices same diag (or triangular same diag) the relations are exactly

   p u s = p v s

   where p and s contains both x and y

2. For band 3x3 matrices same diag

Conjecture:

1. band matrices

    relations are of the form

    p u s = p v s

    with |p| >= dim and |s| >= dim

    (just use the fact that the first occurrences of factor of length d-1 appear
    at the same positions)

2. triangular matrices:

    relations are of the form

    p u s = p v s

    with |p| > dim and |s| > dim
"""

from sage.geometry.polyhedron.parent import Polyhedra
from semigroup_tools import products_n, products_p, is_relation

###################
# Relation finder #
###################

def pretty_relation_string(r1, r2, letters='mn'):
    r"""
    A nice string to represent a relation

    EXAMPLES::

        sage: pretty_relation_string([0,1,1,0],[1,1,0,0], 'xy')
        '(xyy)x = (yyx)x'
        sage: pretty_relation_string([0,0,1,0],[0,1,0,0], 'xy')
        'x(xy)x = x(yx)x'
    """
    p = 0
    while r1[p] == r2[p]: p += 1
    s = -1
    while r1[s] == r2[s]: s -= 1
    s += 1
    return '{}({}){} = {}({}){}'.format(
       ''.join(letters[i] for i in r1[:p]),
       ''.join(letters[i] for i in r1[p:s]),
       ''.join(letters[i] for i in r1[s:]),
       ''.join(letters[i] for i in r2[:p]),
       ''.join(letters[i] for i in r2[p:s]),
       ''.join(letters[i] for i in r2[s:]))

def relations_band(dim, start=3, num_mat=20, filename=None):
    r"""
    List the relations for the band matrices (with identical diagonals) in a
    given dimension.

    INPUT:

    - ``dim`` -- the dimension

    - ``start`` -- the length of product we consider first (default to ``1``)

    - ``num_mat`` -- the number of integer matrices used to check the relation.
      Note that if all the integer matrices satisfy the relations, a (costly)
      symbolic check is performed to guarantee that the relation is satisfied by
      any pair of matrices.

    EXAMPLES::

        sage: relations_tri(2)
        # triangular relations for n = 1
        # triangular relations for n = 2
        # triangular relations for n = 3
        # triangular relations for n = 4
        # triangular relations for n = 5
        # triangular relations for n = 6
        # triangular relations for n = 7
        # triangular relations for n = 8
        # triangular relations for n = 9
        # triangular relations for n = 10
        #  relations (5,5)
        xyyx(xy)xyyx = xyyx(yx)xyyx
        xyyx(xy)yxxy = xyyx(yx)yxxy
        xyyx(yx)xyyx = xyyx(xy)xyyx
        xyyx(yx)yxxy = xyyx(xy)yxxy
        # triangular relations for n = 11
        #  relations (4,7)
        xyyyyx(xy)yxy = xyyyyx(yx)yxy
        xyyyyx(yx)yxy = xyyyyx(xy)yxy
        ...
    """
    # the symbolic max-plus matrices
    ab = symbolic_max_plus_matrices_band(dim, 2)
    one = symbolic_max_plus_identity(dim, ab[0].num_variables())

    # the integer max-plus matrices
    pairs = [random_integer_max_plus_matrices_band(dim, -2**30, 2**30) for _ in range(num_mat)]
    A,B = pairs[0]

    one_int = integer_max_plus_matrix_identity(dim)

    # n     : length of the product
    # i1,m1AA,m1AB : data for the first word
    # i2,m2AA,m2AB : data for the second word
    n = start
    while True:
        if filename is None:
            from sys import stdout
            f = stdout
        else:
            f = open(filename.format(dim,n), 'w')

        # we restrict to products that start and ends with the same letters
        # up to symmetry we can restrict to
        #   X ... X = X ... X
        #   X ... Y = X ... Y

        f.write("# band similar diagonal relations for n = {}\n".format(n))
        print "n = {}".format(n)
        for i1,m1 in products_n(A, B, n-2, one_int):
            m1AA = A * m1 * A
            m1AB = A * m1 * B

            # TODO: fill the things that we know
            # the first and last occurrences of any subword of length d-1 appear
            # in the same position in the left and right members
            for i2,m2 in products_n(A, B, n-2, one_int):
                m2AA = A * m2 * A
                m2AB = A * m2 * B

                if i1 == i2:
                    break

                # we have a combinatorial condition to check (comes from
                # dimension 2):
                # the longest common prefix and longest common suffix must
                # contain the two letters
                j = 0
                while i1[j] == i2[j] == 0:
                    j += 1
                if i1[j] != i2[j]:
                    continue

                # here we first test equality between m1 and m2
                # then we test the relations on all matrices in mats
                # then we check formally using symbolic matrices
                ii1 = [0] + i1 + [0]
                ii2 = [0] + i2 + [0]
                if m1AA == m2AA and \
                   is_relation(ii1, ii2, pairs) and \
                   prod(ab[x] for x in ii1) == prod(ab[x] for x in ii2):
                       f.write(pretty_relation_string(ii1,ii2,'xy'))
                       f.write('\n')
                       f.flush()
                ii1 = [0] + i1 + [1]
                ii2 = [0] + i2 + [1]
                if m1AB == m2AB and \
                   is_relation(ii1, ii2, pairs) and \
                   prod(ab[x] for x in ii1) == prod(ab[x] for x in ii2):
                       f.write(pretty_relation_string(ii1,ii2,'xy'))
                       f.write('\n')
                       f.flush()

        if filename is not None:
            f.close()
        n += 1

def relations_tri_sim_diag(dim, start=3, num_mat=20, filename=None):
    r"""
    List the relations for the upper triangular matrices in a given dimension.

    INPUT:

    - ``dim`` -- the dimension

    - ``start`` -- the length of product we consider first (default to ``1``)

    - ``num_mat`` -- the number of integer matrices used to check the relation.
      Note that if all the integer matrices satisfy the relations, a (costly)
      symbolic check is performed to guarantee that the relation is satisfied by
      any pair of matrices.

    EXAMPLES::

        sage: relations_tri(2)
        # triangular relations for n = 1
        # triangular relations for n = 2
        # triangular relations for n = 3
        # triangular relations for n = 4
        # triangular relations for n = 5
        # triangular relations for n = 6
        # triangular relations for n = 7
        # triangular relations for n = 8
        # triangular relations for n = 9
        # triangular relations for n = 10
        #  relations (5,5)
        xyyx(xy)xyyx = xyyx(yx)xyyx
        xyyx(xy)yxxy = xyyx(yx)yxxy
        xyyx(yx)xyyx = xyyx(xy)xyyx
        xyyx(yx)yxxy = xyyx(xy)yxxy
        # triangular relations for n = 11
        #  relations (4,7)
        xyyyyx(xy)yxy = xyyyyx(yx)yxy
        xyyyyx(yx)yxy = xyyyyx(xy)yxy
        ...
    """
    # the symbolic max-plus matrices
    ab = symbolic_max_plus_matrices_tri_sim_diag(dim, 2)
    one = symbolic_max_plus_identity(dim, ab[0].num_variables())

    # the integer max-plus matrices
    pairs = [random_integer_max_plus_matrices_tri_sim_diag(dim, -2**30, 2**30) for _ in range(num_mat)]
    A,B = pairs[0]

    one_int = integer_max_plus_matrix_identity(dim)

    # n     : length of the product
    # i1,m1AA,m1AB : data for the first word
    # i2,m2AA,m2AB : data for the second word
    n = start
    while True:
        if filename is None:
            from sys import stdout
            f = stdout
        else:
            f = open(filename.format(dim,n), 'w')

        # we restrict to products that start and ends with the same letters
        # up to symmetry we can restrict to
        #   X ... X = X ... X
        #   X ... Y = X ... Y

        f.write("# triangular similar diagonal relations for n = {}\n".format(n))
        print "n = {}".format(n)
        relations = []
        for i1,m1 in products_n(A, B, n-2, one_int):
            m1AA = A * m1 * A
            m1AB = A * m1 * B

            for i2,m2 in products_n(A, B, n-2, one_int):
                m2AA = A * m2 * A
                m2AB = A * m2 * B

                if i1 == i2:
                    break

                # we have a combinatorial condition to check (comes from
                # dimension 2):
                # the longest common prefix and longest common suffix must
                # contain the two letters
                j = 0
                while i1[j] == i2[j] == 0:
                    j += 1
                if i1[j] != i2[j]:
                    continue

                # here we first test equality between m1 and m2
                # then we test the relations on all matrices in mats
                # then we check formally using symbolic matrices
                ii1 = [0] + i1 + [0]
                ii2 = [0] + i2 + [0]
                if m1AA == m2AA and \
                   is_relation(ii1, ii2, pairs) and \
                   prod(ab[x] for x in ii1) == prod(ab[x] for x in ii2):
                       relations.append(pretty_relation_string(ii1,ii2,'xy'))
                ii1 = [0] + i1 + [1]
                ii2 = [0] + i2 + [1]
                if m1AB == m2AB and \
                   is_relation(ii1, ii2, pairs) and \
                   prod(ab[x] for x in ii1) == prod(ab[x] for x in ii2):
                       relations.append(pretty_relation_string(ii1,ii2,'xy'))

        for r in relations:
            f.write(r)
            f.write('\n')
        f.flush()

        if filename is not None:
            f.close()
        n += 1

def relations_tri(dim, start=1, num_mat=10, filename=None):
    r"""
    List the relations for the upper triangular matrices in a given dimension.

    INPUT:

    - ``dim`` -- the dimension

    - ``start`` -- the length of product we consider first (default to ``1``)

    - ``num_mat`` -- the number of integer matrices used to check the relation.
      Note that if all the integer matrices satisfy the relations, a (costly)
      symbolic check is performed to guarantee that the relation is satisfied by
      any pair of matrices.

    EXAMPLES::

        sage: relations_tri(2)
        # triangular relations for n = 1
        # triangular relations for n = 2
        # triangular relations for n = 3
        # triangular relations for n = 4
        # triangular relations for n = 5
        # triangular relations for n = 6
        # triangular relations for n = 7
        # triangular relations for n = 8
        # triangular relations for n = 9
        # triangular relations for n = 10
        #  relations (5,5)
        xyyx(xy)xyyx = xyyx(yx)xyyx
        xyyx(xy)yxxy = xyyx(yx)yxxy
        xyyx(yx)xyyx = xyyx(xy)xyyx
        xyyx(yx)yxxy = xyyx(xy)yxxy
        # triangular relations for n = 11
        #  relations (4,7)
        xyyyyx(xy)yxy = xyyyyx(yx)yxy
        xyyyyx(yx)yxy = xyyyyx(xy)yxy
        ...
    """
    # the symbolic max-plus matrices
    ab = symbolic_max_plus_matrices_tri(dim, 2)
    one = symbolic_max_plus_identity(dim, ab[0].num_variables())

    # the integer max-plus matrices
    mats = [random_integer_max_plus_matrix_tri(dim, -50*dim*dim, 50*dim*dim) for _ in range(num_mat)]
    pairs = [(m1,m2) for m1 in mats for m2 in mats if m1 is not m2]
    one_int = integer_max_plus_matrix_identity(dim)

    # n     : length of the product
    # i1,m1 : data for the first word
    # i2,m2 : data for the second word
    n = start
    while True:
        if filename is None:
            from sys import stdout
            f = stdout
        else:
            f = open(filename.format(dim,n), 'w')

        # TODO (see also in band relations): we can look only at the pair of
        # words
        f.write("# triangular relations for n = {}\n".format(n))
        for k in range(1,n):
            relations = []
            for i1,m1 in products_p(mats[0], mats[1], k-1, n-k, one_int):
                m1 = mats[0] * m1
                i1 = [0] + i1

                for i2,m2 in products_p(mats[0], mats[1], k-1, n-k, one_int):
                    m2 = mats[0] * m2
                    i2 = [0] + i2

                    if i1 == i2:
                        break

                    # here we first test equality between m1 and m2
                    # then we test the relations on all matrices in mats
                    # then we check formally using symbolic matrices
                    if m1 == m2 and \
                       is_relation(i1, i2, pairs) and \
                       prod(ab[x] for x in i1) == prod(ab[x] for x in i2):
                           relations.append(pretty_relation_string(i1,i2,'ab'))

#                for i2,m2 in products(mats[0], mats[1], k, n-k-1, one_int):
#                    m2 = mats[1] * m2
#                    i2 = [1] + i2
#
#                    # here we first test equality between m1 and m2
#                    # then we test the relations on all matrices in mats
#                    # then we check formally using symbolic matrices
#                    if m1 == m2 and \
#                       is_relation(i1, i2, pairs) and \
#                       prod(ab[x] for x in i1) == prod(ab[x] for x in i2):
#                           relations.append(pretty_relation_string(i1,i2))

            if relations:
                f.write("#  relations ({},{})\n".format(k,n-k))
                for r in relations:
                    f.write(r)
                    f.write('\n')
                f.flush()

        if filename is not None:
            f.close()
        n += 1

def filter_relations(r, mats, ab, one):
    ans = []
    for r1,r2 in r:
        if is_relation(r1, r2, mats) and \
            prod(ab[x] for x in r1) == prod(ab[x] for x in r2):
            ans.append(pretty_relation_string(r1,r2))
    return ans

#####
#
####

def brute_force_fibo(dim):
    dim = int(dim)
    w = list(words.FibonacciWord([0,1])[:1000])
    n = 1

    mats = random_integer_max_plus_matrices_tri_sim_diag(dim, -2**20, 2**20)

    while True:
        u = prod(mats[w[i]] for i in range(n))
        v = prod(mats[w[i]] for i in range(n-1,-1,-1))

        if v*mats[0]*mats[1]*u != v*mats[1]*mats[0]*u:
            n += 1
            print n
        else:
            mats = random_integer_max_plus_matrices_tri_sim_diag(dim, -2**20, 2**20)



##############################
# Symbolic max-plus matrices #
##############################

def symbolic_max_plus_identity(dim, N):
    r"""
    Return the ``dim x dim`` identity matrices on ``N`` variables.
    """
    P = Polyhedra(ZZ, N)
    e = P([[],[],[]], None)
    zero = P()

    data = [[zero if i == j else e for j in range(dim)] for i in range(dim)]
    return BinaryMaxPlusMatrix(dim, data)

def symbolic_max_plus_matrices(dim, nb):
    r"""
    Return ``nb`` independent symbolic matrices in dimension ``dim``
    """
    from sage.geometry.polyhedron.parent import Polyhedra

    N = nb * dim * dim
    P = Polyhedra(ZZ, N)
    e = P([[],[],[]], None)
    o = 0
    z = [0]*N
    matrices = []
    for i in range(nb):
        mat = []
        for j in range(dim):
            row = []
            for k in range(dim):
                z[o] = 1
                o += 1
                row.append(P([[z],[],[]],None))
                z[o-1] = 0
            mat.append(row)
        matrices.append(BinaryMaxPlusMatrix(dim, mat))
    return matrices

def symbolic_max_plus_matrices_tri(dim, nb):
    r"""
    Return ``nb`` independent triangular matrices in dimension ``dim``.

    INPUT:

    - ``dim`` -- dimension

    - ``nb`` -- number of matrices
    """
    from sage.geometry.polyhedron.parent import Polyhedra

    N = nb * dim*(dim+1)//2
    P = Polyhedra(ZZ, N)
    e = P([[],[],[]], None)
    o = 0
    z = [0]*N
    matrices = []
    for i in range(nb):
        mat = []
        for j in range(dim):
            row = [e] * j
            for k in range(dim-j):
                z[o] = 1
                o += 1
                row.append(P([[z], [], []], None))
                z[o-1] = 0
            mat.append(row)
        matrices.append(BinaryMaxPlusMatrix(dim, mat))
    return matrices

def symbolic_max_plus_matrices_tri_sim_diag(dim, nb):
    r"""
    Return a set of ``nb`` symbolic matrices with the same diagonal.

    INPUT:

    - ``dim`` -- the dimension

    - ``nb`` -- the number of matrices
    """
    N = nb * dim*(dim-1)//2 + dim
    P = Polyhedra(ZZ, N)
    e = P([[],[],[]], None)

    o = 0
    z = [0]*N
    diag = []
    for k in range(dim):
        z[o] = 1
        o += 1
        diag.append(P([[z], [], []], None))
        z[o-1] = 0

    matrices = []
    for i in range(nb):
        mat = []
        for j in range(dim):
            row = [e] * j
            row.append(diag[j])
            for k in range(dim-j-1):
                z[o] = 1
                o += 1
                row.append(P([[z], [], []], None))
                z[o-1] = 0
            mat.append(row)
        matrices.append(BinaryMaxPlusMatrix(dim, mat))
    return matrices

def symbolic_max_plus_matrices_band(dim, nb):
    r"""
    Return a set of ``nb`` symbolic matrices with the same diagonal.

    INPUT:

    - ``dim`` -- the dimension

    - ``nb`` -- the number of matrices
    """
    N = dim + nb * (dim - 1)
    P = Polyhedra(ZZ, N)
    e = P([[],[],[]], None)

    o = 0
    z = [0]*N
    diag = []
    for k in range(dim):
        z[o] = 1
        o += 1
        diag.append(P([[z], [], []], None))
        z[o-1] = 0

    matrices = []
    for i in range(nb):
        mat = []
        for j in range(dim):
            row = [e] * j
            row.append(diag[j])
            if j < dim-1:
                z[o] = 1
                o += 1
                row.append(P([[z], [], []], None))
                z[o-1] = 0
                row.extend([e] * (dim-j-2))
            mat.append(row)
        matrices.append(BinaryMaxPlusMatrix(dim, mat))
    return matrices



class BinaryMaxPlusMatrix(SageObject):
    r"""
    A symbolic max plus matrix.
    """
    def __init__(self, n, data):
        self._n = n
        if len(data) != n or any(len(x) != n for x in data):
            raise ValueError
        self._data = tuple(tuple(x) for x in data)

    def num_variables(self):
        r"""
        Return the number of variables used in this matrix.
        """
        return self._data[0][0].ambient_space().dimension()

    def equal_coefficients(self, other):
        r"""
        Return the list of equal coefficients between self and other.
        """
        n = self._n
        return [(i,j) for i in range(n) for j in range(n) \
                if self._data[i][j] == other._data[i][j]]

    def __eq__(self, other):
        r"""
        Equality test
        """
        return (self._n == other._n and
                all(self._data[i][j] == other._data[i][j] \
                        for i in range(self._n) for j in range(self._n)))

    def __ne__(self, other):
        r"""
        Difference test
        """
        return not self.__eq__(other)

    def __mul__(self, other):
        r"""
        Multiplication of matrices
        """
        n = self._n
        assert self._n == other._n
        new_data = [[None]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                vertices = []
                for k in range(n):
                    vertices.extend([xi+yi for xi,yi in zip(x,y)] \
                            for x in self._data[i][k].vertices_list() \
                            for y in other._data[k][j].vertices_list())
                new_data[i][j] = self._data[0][0].parent()([vertices,[],[]],None)

        return BinaryMaxPlusMatrix(self._n, new_data)

    def _repr_(self):
        r"""
        String when the object is printed
        """
        return 'A {}x{} max plus matrix'.format(self._n,self._n)

    def __getitem__(self, data):
        r"""
        To obtain an item of the matrix
        """
        i,j = data
        return self._data[i][j].vertices_list()

