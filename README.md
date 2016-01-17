# Sage library for max-plus matrix identities

## Introduction

To run the library you need to install Sage (http://sagemath.org). Then clone
this repository (or download the files `max_plus.py` and `int_max_plus.pyx`).
Go to the directory where you cloned or downloaded the file. Start Sage and run

    sage: %runfile max_plus.py

## Symbolic and integer max plus matrices

Once you are able to load the file `max_plus.py`, you can create symbolic
matrices with the following functions

    symbolic_max_plus_matrices_band(d, n, diag, surdiag, ch, sym)
	symbolic_max_plus_matrices_upper(d, n, diag, surdiag, ch, sym)
    symbolic_max_plus_matrices(d, n, ch, sym)

for respectively band, upper triangular and full matrices. The arguments are

- `d` - the dimension

- `n` - the number of matrices

- `diag` - describes whether the diagonal is zero (`'z'`), constant (`'c'`),
  the same in each of the `n` matrices (`'s'`) or variable (`'v'`). This is
  optional and default to `'v'`.

- `surdiag` - describes whether the surdiagonal is constant, zero, same or
  variable.  This is also optional.

- `ch` - (optional) set the convex hull engine (could be one of 'ppl', 'cdd'
  or 'PALP'). It defaults to 'ppl' which seems to be the fastest.

- `sym` - (optional, default to `True`) specifies the implementation: if `sym`
  is `True` the matrices do not store all the coefficients (only two in the case
  of full matrices and `d` in the case of triangular ones). The multiplication
  in these case should be a little bit faster (by a factor `d^2` for dense
  matrix but much less for triangular ones). This should be set to `False` only
  to debug the code or to play with `diag=c` or `surdiag=c` that are not
  compatible with this option.

For example you can do

    sage: x,y = symbolic_max_plus_matrices_band(3, 2, 'z', 'v')
    sage: x
    A 3x3 symbolic max plus matrix on 4 variables
    sage: print x
	[   0  x0 -oo ]
	[ -oo   0  x1 ]
	[ -oo -oo   0 ]
	sage: y
    A 3x3 symbolic max plus matrix on 4 variables
	sage: print y
	[   0  x2 -oo ]
	[ -oo   0  x3 ]
	[ -oo -oo   0 ]
    sage: p = x*y*x
    A 3x3 symbolic max plus matrix on 4 variables
    sage: print p
	[   0 max(x2, x0) max(x1+x2, x0+x3, x0+x1) ]
	[ -oo           0              max(x3, x1) ]
	[ -oo         -oo                        0 ]
	sage: x*y*x*y == y*x*y*x == x*y*y*x == y*x*x*y
	True

	sage: x2,y2,z2 = symbolic_max_plus_matrices_upper(3, 3)
    sage: x2.num_vars()
    18
	sage: print x2
	[  x0   x3  x5 ]
	[ -oo   x1  x4 ]
	[ -oo  -oo  x2 ]
	sage: print y2
	[  x6   x9  x11 ]
	[ -oo   x7  x10 ]
	[ -oo  -oo   x8 ]
	sage: print z2
	[ x12  x15  x17 ]
	[ -oo  x13  x16 ]
	[ -oo  -oo  x14 ]

And with full 3x3 matrices

    sage: x,y,z = symbolic_max_plus_matrices(3, 3)
    sage: print x
	[ x0 x1 x2 ]
	[ x3 x4 x5 ]
	[ x6 x7 x8 ]
    sage: print y
	[  x9 x10 x11 ]
	[ x12 x13 x14 ]
	[ x15 x16 x17 ]
	sage: print x*x
	[   max(x2+x6, x1+x3, 2x0) max(x2+x7, x1+x4, x0+x1) max(x2+x8, x1+x5, x0+x2) ]
	[ max(x5+x6, x3+x4, x0+x3)   max(x5+x7, 2x4, x1+x3) max(x5+x8, x4+x5, x2+x3) ]
	[ max(x6+x8, x3+x7, x0+x6) max(x7+x8, x4+x7, x1+x6)   max(2x8, x5+x7, x2+x6) ]

From a symbolic matrix you can obtain an integer max plus matrix using the
method `eval`. For that you first need to compile `int_max_plus.pyx`

    sage: %runfile int_max_plus.pyx
    Compiling ./int_max_plus.pyx...
    sage: x,y = symbolic_max_plus_matrices_band(3, 2, 'z', 'v')
    sage: xv = x.eval((1,-1,0,3))
    sage: xv
    [   0   1 -oo ]
	[ -oo   0  -1 ]
	[ -oo -oo   0 ]
    sage: yv = y.eval((1,-1,0,3))
    sage: yv
    [   0   0 -oo ]
	[ -oo   0   3 ]
	[ -oo -oo   0 ]
    sage: p = x*y*x
    sage: pv = p.eval((1,-1,0,3))
    sage: pv
	[   0   1 4 ]
	[ -oo   0 3 ]
	[ -oo -oo 0 ]
    sage: pv == xv*yv*xv
    True

Computations with integer matrices are infinitely faster. Hence, if you intend
to test some relations it is adviced to first test them on a sample of integer
matrices.

## Combinatorics

There are currently few functions to check some relations combinatorially in `B^{sv}_d`.

- `def occurrences(w, u)`: return the occurrences of `u` in `w`

- `def extremal_occurrences(w, u)`: return the 1-extremal occurrences of `u` in
  `w` (the convex hull of this set is the same as all occurrences)

- `def extremal_occurrences2(w,u)`: (experimental) return the extremal
  occurrences of `u` in `w`

- `def extremal_mid_occurrences(p, s, u)`: return the occurrences of `u` in
  `p*s` where the occurrence uses the joker letter `*`.

To check identities you can use the following functions (note that identities
must be written with the letters 'x' and 'y'):

- `is_sv_identity(p, s, d)`: check whether `(pxs,pys)` is an identity in `B^{sv}_d`.

- `is_sv_identity_parallel(p, s, d, prefix_length)`: the same as above but
  with parallelization. The argument `prefix_length` is used to chunk the
  subwords into different jobs: if it is set to `k` then there will be
  `2^k` jobs which correspond to the `2^k` possible prefixes.

Here are some examples

    sage: p = 'xyxxyy'
    sage: s = 'xxyyxy'
    sage: is_sv_identity(p, s, 3)
    True

    sage: p = 'xyyxxxyyyyxxxx'
    sage: s = 'yyyyxxxxyyyxxy'
    sage: %time is_sv_identity(p, s, 5)
	CPU times: user 136 ms, sys: 4 ms, total: 140 ms
	Wall time: 130 ms
	True

	sage: p = 'xyyxxxyyyyxxxxxyyyyy'
	sage: s = 'xxxxxyyyyyxxxxyyyxxy'
	sage: %time is_sv_identity(p, s, 6)
	CPU times: user 812 ms, sys: 8 ms, total: 820 ms
	Wall time: 786 ms
	True
    sage: is_sv_identity(p, s, 7)
    False

	sage: p,s = vincent_sv_prefix_suffix(7)
	sage: is_sv_identity_parallel(p, s, 7, 3, verbose=True)
	PoolWorker-17: new job at 16:43:7
	  ('xyyxxxyyyyxxxxxyyyyyyxxxxxx', 'yyyyyyxxxxxxyyyyyxxxxyyyxxy', 7, ('x', 'x', 'x'))
	PoolWorker-18: new job at 16:43:7
	  ('xyyxxxyyyyxxxxxyyyyyyxxxxxx', 'yyyyyyxxxxxxyyyyyxxxxyyyxxy', 7, ('x', 'x', 'y'))
	PoolWorker-19: new job at 16:43:7
	  ('xyyxxxyyyyxxxxxyyyyyyxxxxxx', 'yyyyyyxxxxxxyyyyyxxxxyyyxxy', 7, ('x', 'y', 'x'))
	PoolWorker-20: new job at 16:43:8
	  ('xyyxxxyyyyxxxxxyyyyyyxxxxxx', 'yyyyyyxxxxxxyyyyyxxxxyyyxxy', 7, ('x', 'y', 'y'))
	PoolWorker-19: job done in 1.26115489006 seconds
	PoolWorker-19: new job at 16:43:9
	  ('xyyxxxyyyyxxxxxyyyyyyxxxxxx', 'yyyyyyxxxxxxyyyyyxxxxyyyxxy', 7, ('y', 'x', 'x'))
	PoolWorker-20: job done in 1.93957805634 seconds
	PoolWorker-20: new job at 16:43:9
	  ('xyyxxxyyyyxxxxxyyyyyyxxxxxx', 'yyyyyyxxxxxxyyyyyxxxxyyyxxy', 7, ('y', 'x', 'y'))
	PoolWorker-18: job done in 2.48127698898 seconds
	PoolWorker-18: new job at 16:43:10
	  ('xyyxxxyyyyxxxxxyyyyyyxxxxxx', 'yyyyyyxxxxxxyyyyyxxxxyyyxxy', 7, ('y', 'y', 'x'))
	PoolWorker-17: job done in 2.97481107712 seconds
	PoolWorker-17: new job at 16:43:10
	  ('xyyxxxyyyyxxxxxyyyyyyxxxxxx', 'yyyyyyxxxxxxyyyyyxxxxyyyxxy', 7, ('y', 'y', 'y'))
	PoolWorker-20: job done in 1.38204503059 seconds
	PoolWorker-19: job done in 2.0873029232 seconds
	PoolWorker-18: job done in 1.68425393105 seconds
	PoolWorker-17: job done in 1.7091588974 seconds
	computation with 4 cpus performed in 4.71341395378 seconds
	True

## Experimental

There is an implementation of barycentric decomposition with PPL Linear
programming. In order to use it you need to install the [pplpy
package](https://pypi.python.org/pypi/pplpy/). It is essentially a fork of some
files in Sage with several bug fixes and improvements. To install it just do

    $ sage -pip install pplpy

Once this is done, you can use the following function

- `barycentric_coordinates(pts, q)`: return the barycentric coordinates of the point `q` with
  respect to the points `pts`. If `q` does not belong to the convex hull of `pts` the value
  `None` is returned instead

(if you did not installed the package, running this function will crash Sage!)

For example

    sage: barycentric_coordinates([(2,3,6),(2,3,18),(15,16,18)], (8,9,13))
	[5/12, 19/156, 6/13]

## Contact

Any comment or remark is welcome at vincentDOTdelecroixATlabriDOTfr
