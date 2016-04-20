from sage_import import *

import itertools
from combinat import extremal_occurrences
from convex_hull import ppl_polytope

def is_vv_identity(left, right, d, prefix=(), status=False):
    r"""
    EXAMPLES::

        sage: from max_plus import *

        sage: p,s = vincent_sv_prefix_suffix(4)
        sage: u = (p+'x'+s).replace('x', 'ab').replace('y', 'ba')
        sage: v = (p+'y'+s).replace('x', 'ab').replace('y', 'ba')
        sage: is_vv_identity(u, v, 4)
        True

        sage: u = (p+'x'+s).replace('x','a').replace('y','b')
        sage: v = (p+'y'+s).replace('x','a').replace('y','b')
        sage: is_vv_identity(u, v, 4)
        False
    """
    n = len(left)
    if n != len(right):
        raise RuntimeError("the two sides must have same length")

    if status:
        output = ''

    # now iterate through all words `u` with given prefix
    pref = tuple(prefix)
    for q in itertools.product('ab', repeat=d-1-len(prefix)):
        u = prefix+q

        # compute the polytope of occurrences and check that it contains the
        # "middle occurrences"
        # TODO: count the number of points in extremal_mid_occurrences and add
        # it to status
        # similarly, if there is a bad descent word, add it
        left_occ = set(tuple(left[:i].count('a') for i in occ) + \
                       tuple(left[:i].count('b') for i in occ) \
                       for occ in extremal_occurrences(left, u))
        right_occ = set(tuple(right[:i].count('a') for i in occ) + \
                       tuple(right[:i].count('b') for i in occ) \
                       for occ in extremal_occurrences(right, u))

        inter = left_occ.intersection(right_occ)
        union = left_occ.difference(right_occ)

        if status:
            t0 = time()
        P = ppl_polytope(inter)
        if status:
            output += 'u = {}\n'.format(''.join(u))
            output += 'num ext. occ.: {}\n'.format(len(inter))
            output += 'num int. occ.: {}\n'.format(len(union))
            output += 'num faces    : {}\n'.format(len(P.minimized_constraints()))
            output += 'num verts    : {}\n'.format(len(P.minimized_generators()))
            output += 'polytope computation in {}secs\n'.format(time()-t0)
            t0 = time()
        for o in union:
            pt = C_Polyhedron(point(Linear_Expression(o,0)))
            if not P.contains(pt):
                if status:
                    output += 'bad occurrence {}\n'.format(o, ''.join(left[i] for i in o))
                return (False,output) if status else False
        if status:
            output += 'containance test in {}secs\n'.format(time()-t0)
            output += '\n'
    return (True,output) if status else True

