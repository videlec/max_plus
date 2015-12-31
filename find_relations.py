from semigroup_tools import (products_n, products_p,
        constraints_from_subwords, minimal_all_subwords_prefix,
        minimal_all_subwords_suffix)

######################
# Experimental stuff #
######################
def prod_symbolic(w, mats):
    r"""
    Recursive computation that tries to minize the number of products

    DO NOT USE... IT DOES NOT SEEMS FASTER!
    """
    assert max(w) <= len(mats)
    if len(w) <= 3 or 3*len(mats) >= 2*len(w):
        return prod(mats[i] for i in w)

    ww = [(w[i],w[i+1]) for i in range(0,len(w)-1,2)]
    P2 = list(set(ww))
    P2_index = {p:i for i,p in enumerate(P2)}
    mmats = [mats[p[0]] * mats[p[1]] for p in P2]

    ww = [P2_index[u] for u in ww]
    if len(w)%2:
        mmats.append(mats[w[-1]])
        ww.append(len(mmats)-1)
    return prod_symbolic(ww, mmats)

def relations_band(dim, start=3, num_mat=5000, limit_coeff=1073741824, filename=None):
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

    .. NOTE::

        All the time in spent in symbolic verification. The only way to make it
        fast is to have a huge ``num_mat`` in order to eliminate the maximum of
        relations without symbolic computations.

    EXAMPLES::

        sage: relations_band(3)
        # band similar diagonal relations for n = 11
        n = 11
        xxyyx(y)xxyxy = xxyyx(x)xxyxy
        xxyyx(y)xxyyx = xxyyx(x)xxyyx
        xxyyx(y)xyyxx = xxyyx(x)xyyxx
        xxyyx(y)yxxyy = xxyyx(x)yxxyy
        xxyyx(y)yyxxy = xxyyx(x)yyxxy
        xyxyy(y)yxxyy = xyxyy(x)yxxyy
        xyxyy(y)yyxxy = xyxyy(x)yyxxy
        xyxyy(y)yyxyx = xyxyy(x)yyxyx
        xyyxx(y)xxyxy = xyyxx(x)xxyxy
        xyyxx(y)xxyyx = xyyxx(x)xxyyx
        xyyxx(y)xyyxx = xyyxx(x)xyyxx
        xyyxx(y)yxxyy = xyyxx(x)yxxyy
        xyyxx(y)yyxxy = xyyxx(x)yyxxy
        # band similar diagonal relations for n = 12
        n = 12
        xxxyyx(y)xxyxy = xxxyyx(x)xxyxy
        xxxyyx(y)xxyyx = xxxyyx(x)xxyyx
        xxxyyx(y)xyyxx = xxxyyx(x)xyyxx
        xxxyyx(y)yxxyy = xxxyyx(x)yxxyy
        xxxyyx(y)yyxxy = xxxyyx(x)yyxxy
        xxyxyx(y)xxyxy = xxyxyx(x)xxyxy
        xxyxyx(y)xxyyx = xxyxyx(x)xxyyx
        xxyxyx(y)xyyxx = xxyxyx(x)xyyxx
        xxyxyy(y)xxyyx = xxyxyy(x)xxyyx
        xxyxyy(y)xyyxx = xxyxyy(x)xyyxx
        xxyxyy(y)yxxyy = xxyxyy(x)yxxyy
        xxyxyy(y)yyxxy = xxyxyy(x)yyxxy
        xxyxyy(y)yyxyx = xxyxyy(x)yyxyx
        xxyyxx(y)xxyxy = xxyyxx(x)xxyxy
        xxyyxx(y)xxyyx = xxyyxx(x)xxyyx
        xxyyxx(y)xyyxx = xxyyxx(x)xyyxx
        xxyyxx(y)yxxyy = xxyyxx(x)yxxyy
        xxyyxx(y)yyxxy = xxyyxx(x)yyxxy
        xxyyx(y)xxxyxy = xxyyx(x)xxxyxy
        xxyyx(yx)xxyxy = xxyyx(xy)xxyxy
        xxyyx(y)xxxyyx = xxyyx(x)xxxyyx
        xxyyx(yx)xxyyx = xxyyx(xy)xxyyx
        xxyyx(y)xxyxxy = xxyyx(x)xxyxxy
        xxyyx(y)xxyxyx = xxyyx(x)xxyxyx
        xxyyx(y)xxyxyy = xxyyx(x)xxyxyy
        xxyyx(y)xxyyxx = xxyyx(x)xxyyxx
        xxyyx(y)xxyyxy = xxyyx(x)xxyyxy
        xxyyx(yx)xyyxx = xxyyx(xy)xyyxx
        xxyyx(y)xxyyyx = xxyyx(x)xxyyyx
        xxyyx(y)xyxxyy = xxyyx(x)xyxxyy
        xxyyx(yx)yxxyy = xxyyx(xy)yxxyy
        xxyyx(y)xyxyxx = xxyyx(x)xyxyxx
        xxyyx(y)xyyxxx = xxyyx(x)xyyxxx
        xxyyx(y)xyyxxy = xxyyx(x)xyyxxy
        xxyyx(yx)yyxxy = xxyyx(xy)yyxxy
        xxyyx(y)xyyyxx = xxyyx(x)xyyyxx
        xxyyx(y)yxxxyy = xxyyx(x)yxxxyy
        xxyyx(yy)xxyyx = xxyyx(xx)xxyyx
        xxyyx(y)yxxyyx = xxyyx(x)yxxyyx
        xxyyx(y)yxxyyy = xxyyx(x)yxxyyy
        xxyyxy(y)xxyyx = xxyyxy(x)xxyyx
        xxyyx(yy)xyyxx = xxyyx(xx)xyyxx
        xxyyx(y)yxyyxx = xxyyx(x)yxyyxx
        xxyyxy(y)xyyxx = xxyyxy(x)xyyxx
        xxyyx(y)yyxxxy = xxyyx(x)yyxxxy
        xxyyx(yy)yxxyy = xxyyx(xx)yxxyy
        xxyyx(y)yyxxyx = xxyyx(x)yyxxyx
        xxyyx(y)yyxxyy = xxyyx(x)yyxxyy
        xxyyxy(y)yxxyy = xxyyxy(x)yxxyy
        ...
    """
    # the symbolic max-plus matrices
    ab = symbolic_max_plus_matrices_band(dim, 2, diag='s', surdiag='v')

    # the integer max-plus matrices
    pairs = [random_integer_max_plus_matrices_band(dim, -limit_coeff, limit_coeff) for _ in range(num_mat)]
    one_int = integer_max_plus_matrix_identity(dim)

    # n     : length of the product
    # ii1   : first word
    # ii2   : second word
    n = start
    while True:
        eliminated_from_int = 0
        eliminated_from_symb = 0
        relations = 0

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
        for i1 in product((0,1), repeat=n-2):
            for pref,suff in ((0,),(0,)),((0,),(1,)):
                ii1 = pref + i1 + suff
                assert len(ii1) == n
                p,s = constraints_from_subwords(ii1, dim)
                if p != -1:
                    i1p = ii1[:p]; i1s = ii1[s:]
                    for i2 in product((0,1), repeat=s-p):
                        ii2 = i1p + i2 + i1s
                        assert len(ii2) == n
                        if ii1 == ii2:
                            break
                        if not is_relation(ii1, ii2, pairs):
                            eliminated_from_int += 1
                            continue
                        if prod(ab[x] for x in ii1) != prod(ab[x] for x in ii2):
                            eliminated_from_symb += 1
                            continue

                        f.write(pretty_relation_string(ii1,ii2,'xy'))
                        f.write('\n')
                        f.flush()
                        relations += 1
                        del i2  # makes itertools faster!
            del i1  # makes itertools faster!

        if filename is not None:
            f.close()
        n += 1

        print "  int elimination  : {}".format(eliminated_from_int)
        print "  symb elimination : {}".format(eliminated_from_symb)
        print "  relations        : {}".format(relations)

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
    abc = symbolic_max_plus_matrices_tri_sim_diag(dim, 3)
    one_ab = symbolic_max_plus_identity(dim, ab[0].num_vars())
    one_abc = symbolic_max_plus_identity(dim, abc[0].num_vars())


    # the integer max-plus matrices
    pairs = [random_integer_max_plus_matrices_tri_sim_diag(dim, -2**30, 2**30) for _ in range(num_mat)]
    A,B = AB = pairs[0]

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

        f.write("# triangular similar diagonal relations for n = {}\n".format(n))
        relations = []

        # we first want to consider the prefix/suffix such that we got a
        # relation whatever is in between (of the same length)
        # how do we identify these relations?

        for p in minimal_all_subwords_prefix(n - 2*(dim-1), dim-1):
            m_prefix = prod(AB[i] for i in p)
            for s in minimal_all_subwords_suffix(n - len(p) - 1, dim-1):
                m_suffix = prod(AB[i] for i in s)
                for i1,m1 in products_n(A,B,n-len(p)-len(s)):
                    for i2,m2 in products_n(A,B,n-len(p)-len(s)):
                        if i1 == i2:
                            break
#
#                # here we first test equality between m1 and m2
#                # then we test the relations on all matrices in mats
#                # then we check formally using symbolic matrices
#                ii1 = [0] + i1 + [0]
#                ii2 = [0] + i1[:ppos] + i2 + i1[spos:] + [0]
#                print "{} =?= {}".format(ii1,ii2)
#                if m1AA == m2AA and \
#                   is_relation(ii1, ii2, pairs) and \
#                   prod(ab[x] for x in ii1) == prod(ab[x] for x in ii2):
#                       f.write(pretty_relation_string(ii1,ii2,'xy'))
#                       f.write('\n')
#                       f.flush()
#                ii1 = [0] + i1 + [1]
#                ii2 = [0] + i1[:ppos] + i2 + i1[spos:] + [0]
#                if m1AB == m2AB and \
#                   is_relation(ii1, ii2, pairs) and \
#                   prod(ab[x] for x in ii1) == prod(ab[x] for x in ii2):
#                       f.write(pretty_relation_string(ii1,ii2,'xy'))
#                       f.write('\n')
#                       f.flush()
#
#        if filename is not None:
#            f.close()

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
    one = symbolic_max_plus_identity(dim, ab[0].num_vars())

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

def relations_band_vc(dim, start=1, num_mat=500, filename=None):
    # the symbolic max-plus matrices
    ab = symbolic_max_plus_matrices_tri(dim, 2)
    one = symbolic_max_plus_identity(dim, ab[0].num_vars())

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

        f.write("# Band^cv_{} relations of length n = {}\n".format(dim,n))

        for n0 in range(1,n-1):
            eliminated_from_int = 0
            eliminated_from_symb = 0
            nb_identities = 0

            relations = []
            n1 = n - n0
            for i1,m1 in products_p(pairs[0][0], pairs[0][1], n0, n1, one_int):
                p,s = constraints_from_subwords(i1, dim)
                if p == -1:
                    continue

                mp = prod(pairs[0][k] for k in i1[:p])
                ms = prod(pairs[1][k] for k in i1[s:])

                nn1 = sum(i1[p:s])
                nn0 = s-p-nn1
                for ii2,mm2 in products_p(pairs[0][0], pairs[0][1], nn0, nn1, one_int):
                    i2 = i1[:p] + ii2 + i1[s:]
                    assert len(i1) == len(i2) and sum(i1) == sum(i2)
                    if i1 == i2:
                        break

                    if not is_relation(tuple(i1), tuple(i2), pairs):
                        eliminated_from_int += 1
                    elif prod(ab[x] for x in i1) != prod(ab[x] for x in i2):
                        eliminated_from_symb += 1
                    else:
                        nb_identities += 1
                        relations.append(pretty_relation_string(i1,i2,'xy'))

            print "({},{})".format(n0,n1)
            print "  int elimination  : {}".format(eliminated_from_int)
            print "  symb elimination : {}".format(eliminated_from_symb)
            print "  identities       : {}".format(nb_identities)

            if relations:
                f.write("#  relations ({},{})\n".format(n0,n1))
                for r in relations:
                    f.write(r)
                    f.write('\n')
                    f.flush()

        if filename is not None:
            f.close()
        n += 1

def conj_min_relations_band_vc(dim, start=1, num_mat=500, filename=None):
    r"""
    Look at relations of the form

    M * x*y * M == M * y*x * M

    for matrices in B^vc where M contains the same number of x and y.
    """
    # the symbolic max-plus matrices
    ab = symbolic_max_plus_matrices_band(dim, 2, diag='v', surdiag='c')
    one = symbolic_max_plus_identity(dim, ab[0].num_vars())
    xy = ab[0]*ab[1]
    yx = ab[1]*ab[0]

    # the integer max-plus matrices
    pairs = [random_integer_max_plus_matrices_band(dim, -50*dim*dim, 50*dim*dim,
        ord('v'), ord('c')) for _ in range(num_mat)]
    one_int = integer_max_plus_matrix_identity(dim)

    # n  : length of m
    # i  : word of m
    # i1 : word of mxym
    # i2 : word of myxm
    n = 1
    while True:
        if filename is None:
            from sys import stdout
            f = stdout
        else:
            f = open(filename.format(dim,n), 'w')

        header = "# Band^vc_{} relations with n = {}".format(dim,n)
        f.write(header)
        f.write("\n")
        if filename:
            print header

        eliminated_from_int = 0
        eliminated_from_symb = 0
        nb_identities = 0
        for i in product_p(n-1, n):
            i = (0,) + tuple(i)
            i1 = i + (0,1) + i
            i2 = i + (1,0) + i

            if not is_relation(i1, i2, pairs, upper=True):
                eliminated_from_int += 1
                continue

            print "  potential relation:", pretty_relation_string(i1,i2,'xy')
            sys.stdout.flush()

            t0 = time.clock()
            m = prod(ab[x] for x in i)
            t1 = time.clock()
            print "    m    computed in {} seconds".format(t1-t0)
            print "    nb pts in convex hull: {}".format([len(m[i,j]) for i in range(dim) for j in range(i,dim)])
            t0 = time.clock()
            mxym = m * xy * m
            t1 = time.clock()
            print "    mxym computed in {} seconds".format(t1-t0)
            print "    nb pts in convex hull: {}".format([len(mxym[i,j]) for i in range(dim) for j in range(i,dim)])
            t0 = time.clock()
            myxm = m * yx * m
            t1 = time.clock()
            print "    myxm computed in {} seconds".format(t1-t0)
            print "    nb pts in convex hull: {}".format([len(myxm[i,j]) for i in range(dim) for j in range(i,dim)])
            if mxym != myxm:
                print "    ... not a relation\n"
                sys.stdout.flush()
                eliminated_from_symb += 1
                continue

            print "    NEW RELATION\n"
            nb_identities += 1
            f.write(pretty_relation_string(i1,i2,'xy'))
            f.write('\n')
            f.flush()
            print

        print "  int elimination  : {}".format(eliminated_from_int)
        print "  symb elimination : {}".format(eliminated_from_symb)
        print "  identities       : {}".format(nb_identities)
        print

        if filename is not None:
            f.close()
        if nb_identities:
            return
        n += 1
        if n == 9:
            break

# TODO
#  relations in B^sv?
#  d-relations for B^cv are {(u,v): Subwords_{d-1}(u) = Subwords_{d-1}(v)}

def filter_relations(r, mats, ab, one):
    ans = []
    for r1,r2 in r:
        if is_relation(r1, r2, mats) and \
            prod(ab[x] for x in r1) == prod(ab[x] for x in r2):
            ans.append(pretty_relation_string(r1,r2))
    return ans

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



