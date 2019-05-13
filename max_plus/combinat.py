from __future__ import print_function, division, absolute_import

import itertools

from .sage_import import *

#########################
# Combinatorial methods #
#########################
def runs(w):
    r"""
    EXAMPLES::

    sage: from max_plus.combinat import runs
    sage: runs([0,0,1,1,1,0])
    [2, 3, 1]
    sage: runs([])
    []
    sage: runs([0]*10)
    [10]
    """
    ans = []
    if not w:
        return ans

    n = len(w)
    i = 0
    while i < n:
        a = w[i]
        j = i+1
        while j < n and w[j] == a:
            j += 1
        ans.append(j-i)
        i = j
    return ans

def has_all_subwords(w, r):
    r"""
    EXAMPLES::

        sage: from max_plus.sv_identities import has_all_subwords
        sage: W = FiniteWords([0,1])
        sage: for w in W.iterate_by_length(4):
        ....:     print(w, has_all_subwords(w, 2))
        0000 False
        0001 False
        0010 False
        0011 False
        0100 False
        0101 True
        0110 True
        0111 False
        1000 False
        1001 True
        1010 True
        1011 False
        1100 False
        1101 False
        1110 False
        1111 False

        sage: sum(has_all_subwords(w,3) for w in W.iterate_by_length(4))
        0
        sage: sum(has_all_subwords(w,3) for w in W.iterate_by_length(5))
        0
        sage: sum(has_all_subwords(w,3) for w in W.iterate_by_length(6))
        8
        sage: sum(has_all_subwords(w,3) for w in W.iterate_by_length(7))
        40
        sage: sum(has_all_subwords(w,3) for w in W.iterate_by_length(8))
        128
    """
    n = len(w)
    i = 0
    k = 0
    while i < n and k < r:
        i0 = i
        i += 1
        while i < n and w[i] == w[i0]:
            i += 1
        if i == n:
            return False
        k += 1
        i += 1

    return k == r

def occurrences(w, u, W=None):
    r"""
    Return the set of occurrences of ``u`` in ``w`` (as subword).

    If the word ``w`` contains one or more jokers (i.e. the letter ``'*'``) then
    a pair of lists is returned. The first one is made of occurrences that does
    not pass through joker and the other one is the complement.

    EXAMPLES::

        sage: from max_plus.combinat import occurrences

        sage: occurrences('abbabab', 'abb')
        [(0, 1, 2), (0, 1, 4), (0, 2, 4), (0, 1, 6), (0, 2, 6),
         (0, 4, 6), (3, 4, 6)]
        sage: w = 'abbabaabbaa'; u = 'abaa'
        sage: all(w[i]+w[j]+w[k]+w[l] == u for (i,j,k,l) in occurrences(w,u))
        True

        sage: s,t = occurrences('xyyxx*yyxxy', 'xyx')
        sage: s
        [(0, 1, 3), (0, 2, 3), (0, 1, 4), ..., (3, 7, 9), (4, 7, 9)]
        sage: t
        [(0, 1, 5), (0, 2, 5), (0, 5, 8), ..., (5, 6, 9), (5, 7, 9)]
    """
    # NOTE: this method is efficient because it parses only once the word w. On
    # the other hand there is a clear waste of memory if we just want to go
    # through the occurrences and not make the list of them (which is *not* our
    # case). Though we could eliminate some of the non-extremal occurrences but
    # that would be some work to implement (?)
    pos = [[] for _ in range(len(u))]
    pos2 = [[] for _ in range(len(u))]

    for i,letter in enumerate(w):
        for j in range(len(u)-1,0,-1):
            if letter == '*':
                pos2[j].extend(x + (i,) for x in pos[j-1])
            elif letter == u[j]:
                pos[j].extend(x + (i,) for x in pos[j-1])
                pos2[j].extend(x + (i,) for x in pos2[j-1])
        if letter == '*':
            pos2[0].append((i,))
        elif letter == u[0]:
            pos[0].append((i,))

    return pos[-1] if not pos2[0] else (pos[-1], pos2[-1])

def letter_extremal_occurrences(w, u):
    r"""
    Return the set of letter-extremal occurrences of ``u`` in ``w``.

    An occurrence is letter-extremal if the letters can not move inside the occurrence.
    This is a subset of all occurrences of ``u`` in ``w`` but *much* that
    defines the same convex hull.

    You should actually look at :func:`extremal_occurrences` which is
    even smarter.

    EXAMPLES::

        sage: from max_plus.combinat import (occurrences,
        ....:     letter_extremal_occurrences)

        sage: letter_extremal_occurrences('aabb','ab')
        [(0, 2), (1, 2), (0, 3), (1, 3)]
        sage: letter_extremal_occurrences('aaabbb','ab')
        [(0, 3), (2, 3), (0, 5), (2, 5)]

        sage: p = 'xyyxxxyyy'
        sage: s = 'xxxyyyxxy'
        sage: o1 = occurrences(p+'x'+s, 'xyx')
        sage: o1
        [(0, 1, 3), (0, 2, 3), (0, 1, 4), ...,   (11, 15, 17), (12, 15, 17)]
        sage: len(o1)
        138
        sage: o2 = letter_extremal_occurrences(p+'x'+s, 'xyx')
        sage: o2
        [(0, 1, 3), (0, 2, 3), (5, 6, 9), ...,  (0, 15, 17), (12, 15, 17)]
        sage: len(o2)
        13
        sage: Polyhedron(o1) == Polyhedron(o2)
        True

        sage: from itertools import product
        sage: for w in ('abbaab', 'aabbaabb', 'abbaaabbbaaaabbbaab',
        ....:           'aaaabbbbaaaabbbbaaaa'):
        ....:     for n in (1,2,3,4):
        ....:         for u in product('ab', repeat=n):
        ....:             P1 = Polyhedron(occurrences(w,u))
        ....:             P2 = Polyhedron(letter_extremal_occurrences(w,u))
        ....:             assert P1 == P2
    """
    if len(u) == 0:
        return [()]

    pos = [[] for _ in range(len(u))]

    letters = set(w)
    n = len(w)

    # 1. find next left and next right occurrences of letters
    next_left = [None]*n
    next_right = [None]*n
    last = {letter: -1 for letter in letters}
    for i,letter in enumerate(w):
        next_left[i] = last[letter]
        last[letter] = i
    last = {letter: n for letter in letters}
    for i,letter in enumerate(reversed(w)):
        next_right[n-i-1] = last[letter]
        last[letter] = n-i-1

    # 2. run through w
    for i,letter in enumerate(w):
        for j in range(len(u)-1,0,-1):
            if letter == u[j]:
                for x in pos[j-1]:
                    if next_left[x[-1]] <= x[-2] or next_right[x[-1]] >= i:
                        pos[j].append(x + (i,))
            else:
                k = 0
                while k < len(pos[j-1]):
                    x = pos[j-1][k]
                    if next_left[x[-1]] > x[-2] and next_right[x[-1]] < i:
                        del pos[j-1][k]
                    else:
                        k += 1
        if letter == u[0]:
            pos[0].append((-1,i))

    return [x[1:] for x in pos[-1] if next_left[x[-1]] <= x[-2] or next_right[x[-1]] >= n]

def extremal_occurrences(w, u, W=None, verbose=False):
    r"""
    Return the set of extremal occurrences of ``u`` in ``w``.

    An occurrence is extremal if the blocks can not move inside the occurrence.
    This is a subset of all occurrences of ``u`` in ``w`` but *much* that
    defines the same convex hull.

    EXAMPLES::

        sage: from max_plus.combinat import (occurrences,
        ....:     extremal_occurrences, letter_extremal_occurrences)

        sage: extremal_occurrences('aaaaa', 'aaa')
        [(0, 1, 2), (0, 1, 4), (0, 3, 4), (2, 3, 4)]
        sage: extremal_occurrences('abababa', 'ab')
        [(0, 1), (0, 5), (4, 5)]
        sage: extremal_occurrences('aabaabaabaa', 'ab')
        [(0, 2), (1, 2), (0, 8), (7, 8)]

        sage: extremal_occurrences([0,0,1,0,1], [0,1])
        [(0, 2), (1, 2), (0, 4), (3, 4)]

        sage: W = FiniteWords([0,1])
        sage: extremal_occurrences(W([0,0,1,0,1]), [0,1], W)
        [(0, 2), (1, 2), (0, 4), (3, 4)]

    Note the difference with ``extremal_occurrences``::

        sage: letter_extremal_occurrences('aaaaa', 'aaa')
        [(0, 1, 2), (0, 2, 3), (1, 2, 3), (0, 1, 4), (1, 2, 4), (0, 3, 4), (2, 3, 4)]
        sage: letter_extremal_occurrences('abababa', 'ab')
        [(0, 1), (2, 3), (0, 5), (4, 5)]
        sage: letter_extremal_occurrences('aabaabaabaa', 'ab')
        [(0, 2), (1, 2), (4, 5), (0, 8), (7, 8)]

    Some more tests::

        sage: extremal_occurrences('aaabcdefghaaa', 'aa')
        [(0, 1), (2, 10), (0, 12), (11, 12)]
        sage: extremal_occurrences('bcdef', 'a')
        []
        sage: extremal_occurrences('bacadeaf', 'a')
        [(1,), (6,)]

    Some check for equality of polyhedra::

        sage: w = 'aaabbbaaabbbaaabbb'
        sage: u = 'aaab'
        sage: o0 = occurrences(w,u)
        sage: o1 = letter_extremal_occurrences(w,u)
        sage: o2 = extremal_occurrences(w,u)
        sage: len(o0), len(o1), len(o2)
        (315, 41, 25)
        sage: Polyhedron(o0) == Polyhedron(o1) == Polyhedron(o2)
        True
        sage: u = 'aabaa'
        sage: o0 = occurrences(w,u)
        sage: o1 = letter_extremal_occurrences(w,u)
        sage: o2 = extremal_occurrences(w,u)
        sage: len(o0), len(o1), len(o2)
        (270, 60, 42)
        sage: Polyhedron(o0) == Polyhedron(o1) == Polyhedron(o2)
        True

        sage: from itertools import product
        sage: for w in ('aaaabbbb', 'abbaab', 'aabbaabb', 'abababababababab',
        ....:           'abbaaabbbaaaabbbaab',
        ....:           'aaaabbbbaaaabbbbaaaa'):
        ....:     for n in (1,2,3,4):
        ....:         for u in product('ab', repeat=n):
        ....:             P0 = Polyhedron(occurrences(w,u))
        ....:             P1 = Polyhedron(letter_extremal_occurrences(w,u))
        ....:             P2 = Polyhedron(extremal_occurrences(w,u))
        ....:             assert P0 == P1 == P2
    """
    if not u:
        return [()]

    # 1. Preprocessing of u: we compute the set of its factors inside a suffix
    # trie. To each factor of u corresponds a number in {0, 1, ..., n-1} where n
    # is the number of factors. The used variables are
    #  tf: the transition function (describe what happens when adding a letter)
    #  sl: the suffix link (describe what happens when removing the left letter)
    #  lengths: the length of factors
    if W is None:
        alphabet = sorted(set(u))
        W = FiniteWords(alphabet)
    else:
        alphabet = W.alphabet()

    S = SuffixTrie(W(u, check=False))
    sl = S._suffix_link    # suffix link
    assert sl[0] == -1
    sl[0] = None           # replace -1 by None
    tf = {}                # transition function
    for (r,letter),s in S._transition_function.iteritems():
        tf[(r,letter[0])] = s
    if verbose:
        print("sl = {}".format(sl))
        print("tf = {}".format(tf))
    lengths  = [None]*len(sl)
    lengths[0] = 0      # lengths of the factors
    todo = [0]
    while todo:
        s0 = todo.pop()
        for letter in alphabet:
            t = (s0,letter)
            if t in tf:
                s1 = tf[t]
                lengths[s1] = lengths[s0] + 1
                todo.append(s1)
    if verbose:
        print("lengths = {}".format(lengths))


    # 2. run through w
    state = 0  # maximal suffix of u at the current position
    fact_occ = [-1] * len(sl)  # occurrences of factors of u
    pos = [[] for _ in range(len(u))]  # extremal subword occurrences of u
                                       # are stored as
                                       # [lb, (i0,j0,f0), (i1,j1,f1), ...]
                                       # where:
                                       #   lb: whether the last block is left
                                       #   blocked
                                       #   i0: starting point of a factor
                                       #   j0: endpoint of a factor
                                       #   f0: the index of the factor in the
                                       #       suffix trie
    for iw,letter in enumerate(w):
        if verbose:
            print("iw = {} letter={}".format(iw,letter))
            print("fact_occ = {}".format(fact_occ))
            print("pos")
            for ll in range(len(pos)):
                print("    {}: {}".format(ll,pos[ll]))
        # a. compute the maximum state that can be extended as well as the new
        #    state after having read letter
        old_state = state
        new_state = tf.get((state,letter))
        while state and new_state is None:
            state = sl[state]
            new_state = tf.get((state,letter))
        if new_state is None:
            new_state = 0
        if verbose:
            print("state = {}, new_state = {}".format(state, new_state))

        # b. process subword occurrences
        for length in range(len(u)-1,0,-1):
            if verbose:
                print("  length={}".format(length))
            if letter == u[length]:
                k = 0
                while k < len(pos[length-1]):
                    x = pos[length-1][k]
                    i0, j0, s0 = x[-1]
                    j1 = x[-2][1] if len(x) > 2 else 0
                    lb = x[0]
                    if verbose: print("  x={}  lb={}  j1={}".format(x,lb,j1))
                    if j0 == iw:
                        # case when the last block can be continued
                        xx = x[:-1]
                        ss = tf[(s0,letter)]
                        xx.append((i0, j0+1, ss))
                        pos[length].append(xx)
                        xx[0] = xx[0] or fact_occ[ss] < j1 + lengths[ss]
                        k += 1
                        if verbose: print("  continue last block {}->{}".format(s0,ss))
                    elif lb or fact_occ[s0] == j0:
                        # case when the last block is either blocked on the left
                        # or on the right
                        xx = x[:]
                        ss = tf[(0,letter)]
                        xx.append((iw, iw+1, ss))
                        xx[0] = fact_occ[ss] < j0 + 1
                        pos[length].append(xx)
                        k += 1
                        if verbose: print("  new block {}".format(ss))
                    else:
                        if verbose: print("  delete x")
                        # we remove x
                        del pos[length-1][k]

        if letter == u[0]:
            ss = tf[(0,letter)]
            xx = [fact_occ[ss] == -1, (iw,iw+1,ss)]
            pos[0].append(xx)

        # update the last occurrences of factors of u and switch to the new
        # state
        state = new_state
        s = state
        while s is not None:
            fact_occ[s] = iw+1
            s = sl[s]
        if verbose: print()

    return [sum((tuple(range(i,j)) for (i,j,k) in x[1:]),()) for x in pos[-1] if x[0] or fact_occ[x[-1][2]] == x[-1][1]]

def prefix_suffix_all_subwords(left, r):
    r"""
    Return a list of the same length as left that is filled with the minimal
    prefix and suffix that contain all subwords of length ``r``.

    Return a pair ``(is_trivial, right)``.

    INPUT:

    - ``left`` -- a list that corresponds to a word on a binary alphabet

    - ``r`` -- a non-negative integer

    EXAMPLES::

        sage: from max_plus.combinat import prefix_suffix_all_subwords

        sage: prefix_suffix_all_subwords(list('ab'), 0)
        (False, [None, None])

        sage: prefix_suffix_all_subwords(list('ab'), 1)
        (True, ['a', 'b'])
        sage: prefix_suffix_all_subwords(list('aba'), 1)
        (True, ['a', 'b', 'a'])
        sage: prefix_suffix_all_subwords(list('abbba'), 1)
        (False, ['a', 'b', None, 'b', 'a'])
        sage: prefix_suffix_all_subwords(list('aabbbaa'),1)
        (False, ['a', 'a', 'b', None, 'b', 'a', 'a'])

        sage: prefix_suffix_all_subwords(list('aabbaabbaa'), 2)
        (True, ['a', 'a', 'b', 'b', 'a', 'a', 'b', 'b', 'a', 'a'])
        sage: prefix_suffix_all_subwords(list('abaababababa'), 2)
        (False, ['a', 'b', 'a', 'a', 'b', None, None, None, 'b', 'a', 'b', 'a'])
        sage: prefix_suffix_all_subwords(list('bbaabababa'), 2)
        (False, ['b', 'b', 'a', 'a', 'b', None, 'b', 'a', 'b', 'a'])
        sage: prefix_suffix_all_subwords(list('abbababbab'), 2)
        (False, ['a', 'b', 'b', 'a', None, 'a', 'b', 'b', 'a', 'b'])
        sage: prefix_suffix_all_subwords(list('abbababbaab'), 2)
        (False, ['a', 'b', 'b', 'a', None, None, None, 'b', 'a', 'a', 'b'])
    """
    n = len(left)

    # minimal prefix
    i = 0
    k = 0
    while i < n and k < r:
        i0 = i
        i += 1
        while i < n and left[i] == left[i0]:
            i += 1
        k += 1
        i += 1
    i -= 1

    if i == n-1:
        return True, list(left)

    # minimal suffix
    j = n-1
    k = 0
    while j >= 0 and k < r:
        j0 = j
        j -= 1
        while j >= 0 and left[j] == left[j0]:
            j -= 1
        k += 1
        j -= 1
    j += 1

    if j <= i+1:
        return True, list(left)
    else:
        return False, list(left[:i+1]) + [None] * (n-(i+1)-(n-j)) + list(left[j:])

def iterate_over_holes(u, v, holes, alphabet):
    r"""
    Fill the hole and check that u > v
    and that the pair u[::-1], v[::-1] is smaller than u,v
    """
    for h in itertools.product(alphabet, repeat=len(holes)):
        for j in range(len(holes)):
            v[holes[j]] = h[j]
        tv = tuple(v)
        if u <= tv:
            return
        iu = u[::-1]
        iv = v[::-1]
        if iu <= u and iv <= u and (iu <= tv or iv <= tv):
            yield tv

