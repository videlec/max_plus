import unittest

import sage.all
from sage.misc.prandom import randint
from sage.misc.misc_c import prod
from sage.combinat.words.words import FiniteWords

from max_plus import *
from max_plus.max_plus_int import (random_integer_max_plus_matrices_band,
        is_relation)

W = FiniteWords((0,1))
w_0 = W((0,))
w_1 = W((1,))
w_01 = W((0,1))
w_10 = W((1,0))


class TestIdentities(unittest.TestCase):
    def all_sv_verifs(self, u, v, dim):
        a,b = set(u).union(set(v))
        a1,b1 = symbolic_max_plus_matrices_band(dim, 2, 's', 'v', typ='sym')
        d1 = {a:a1, b:b1}
        a2,b2 = symbolic_max_plus_matrices_band(dim, 2, 's', 'v', typ='full')
        d2 = {a:a2, b:b2}

        mu1 = prod(d1[i] for i in u)
        mu2 = prod(d2[i] for i in u)
        mv1 = prod(d1[i] for i in v)
        mv2 = prod(d2[i] for i in v)

        is_identity = is_sv_identity(u,v,dim)

        self.assertEqual(is_identity,
                         mu1 == mv1,
                         'u = {} v = {} dim = {}'.format(u,v,dim))

        self.assertEqual(mu1 == mv1, mu2 == mv2,
                         'u = {} v = {} dim = {}'.format(u,v,dim))

        if is_identity:
            elts = [random_integer_max_plus_matrices_band(
                       dim, -1000, 1000, ord('s'), ord('v')) \
                    for _ in range(100)]
            zo = {a:0, b:1}
            t0 = tuple(zo[i] for i in u)
            t1 = tuple(zo[i] for i in v)
            self.assertTrue(is_relation(t0,t1, elts, True),
                            'u = {} v = {} dim = {}'.format(u,v,dim))
            self.assertTrue(is_relation(t0,t1,elts,False),
                            'u = {} v = {} dim = {}'.format(u,v,dim))

    def all_vv_verifs(self, u, v, dim):
        a,b = set(u).union(set(v))
        a1,b1 = symbolic_max_plus_matrices_band(dim, 2, 'v', 'v', typ='sym')
        d1 = {a:a1, b:b1}
        a2,b2 = symbolic_max_plus_matrices_band(dim, 2, 'v', 'v', typ='full')
        d2 = {a:a2, b:b2}

        mu1 = prod(d1[i] for i in u)
        mu2 = prod(d2[i] for i in u)
        mv1 = prod(d1[i] for i in v)
        mv2 = prod(d2[i] for i in v)

        is_identity = is_vv_identity(u,v,dim)

        self.assertEqual(is_identity,
                         mu1 == mv1,
                         'u = {} v = {} dim = {}'.format(u,v,dim))
        self.assertEqual(mu1 == mv1, mu2 == mv2,
                         'u = {} v = {} dim = {}'.format(u,v,dim))
        if is_identity:
            elts = [random_integer_max_plus_matrices_band(
                       dim, -1000, 1000, ord('v'), ord('v')) \
                    for _ in range(100)]
            zo = {a:0, b:1}
            t0 = tuple(zo[i] for i in u)
            t1 = tuple(zo[i] for i in v)
            self.assertTrue(is_relation(t0,t1, elts, True),
                            'u = {} v = {} dim = {}'.format(u,v,dim))
            self.assertTrue(is_relation(t0,t1,elts,False),
                            'u = {} v = {} dim = {}'.format(u,v,dim))

    def test_vincent_identities(self):
        for dim in range(2,4):
            p,s = vincent_sv_prefix_suffix(dim)
            self.all_sv_verifs(p+'x'+s, p+'y'+s, dim)
            p = p.replace('x','ab').replace('y','ba')
            s = s.replace('x','ab').replace('y','ba')
            self.all_vv_verifs(p+'ab'+s,p+'ba'+s,dim)

    def test_vv_random(self):
        for _ in range(10):
            p = W.random_element(randint(4,25))
            s = W.random_element(randint(4,25))
            self.all_vv_verifs(p + w_01 + s, p + w_10 + s, randint(2,3))

    def test_sv_random(self):
        for _ in range(10):
            p = W.random_element(randint(4,19))
            s = W.random_element(randint(4,19))
            self.all_sv_verifs(p + w_0 + s, p + w_1 + s, randint(1,4))
            self.all_sv_verifs(p + w_01 + s, p + w_10 + s, randint(1,4))

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestIdentities)
    unittest.TextTestRunner(verbosity=2).run(suite)
