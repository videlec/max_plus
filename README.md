# Sage library for max-plus matrix identities

## Introduction

To run the library you need to install Sage (http://sagemath.org). Once
Sage is installed on your computer, you need to clone this repository and
install it by running

    $ sage -pip install . [--user]

The option `--user` is optional and install the package in your home
directory.

To check that it works just start Sage and run

    sage: import max_plus

If you do not get error any message with the above line, then you are done.

If you want to install a newer version you will need to specify the `--upgrade`
option

    $ sage -pip install . --upgrade [--user]

## Symbolic max plus matrices

The max plus symbolic module define a certain number of functions

    symbolic_max_plus_matrices_band(d, n, diag, surdiag, ch, typ)
	symbolic_max_plus_matrices_upper(d, n, diag, surdiag, ch, typ)
    symbolic_max_plus_matrices(d, n, ch, typ)
    symbolic_max_plus_identity(d, nvar, ch)

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

- `typ` - (optional, default to `sym`) specifies the implementation: if `'sym'`
  the matrices do not store all the coefficients (only two in the case
  of full matrices and `d` in the case of triangular ones). The multiplication
  in these case should be a little bit faster (by a factor `d^2` for dense
  matrix but much less for triangular ones). This should be set to `'full'` only
  to debug the code or to play with `diag=c` or `surdiag=c` that are not
  compatible with this option.

For example you can do

	sage: from max_plus import *

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
method `eval`:: 

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

## Integer max plus matrices

To create integer matrices you can do as follows

	sage: m1 = max_plus.max_plus_int.IntegerMaxPlusMatrix(2, 2, [0,3,1,2])
	sage: m2 = max_plus.max_plus_int.IntegerMaxPlusMatrix(2, 2, [-5,-1,0,2])
	sage: m1
	[ 0 3 ]
	[ 1 2 ]
	sage: m2
	[ -5 -1 ]
	[  0  2 ]
	sage: m1*m2
	[ 3 5 ]
	[ 2 4 ]

Using a program of Dustin Cartwright (available at http://www.math.utk.edu/~cartwright/rank/)
it is possible to compute the Barvinok rank:

    sage: m = IntegerMaxPlusMatrix(3, 3, [2,3,-2,1,3,5,4,5,1])
    sage: m.barvinok_rank()
    2

Note that it is also possible to set some entries to -infinity

	sage: moo = max_plus.max_plus_int.minus_infinity()
	sage: m3 = max_plus.max_plus_int.IntegerMaxPlusMatrix(2, 2, [0,moo,1,moo])
	sage: m3
	[ 0 -oo ]
	[ 1 -oo ]

The integer matrices are currently immutable (it is not possible to modify them).

## Combinatorics

To check identities you can use the following functions

    is_sv_identity(left, right, d, W, prefix, check_common_factors)
    is_vv_identity(left, right, d, W, prefix, status)
    is_sv_identity_parallel(left, right, d, W, prefix_length, ncpus, verbose, check_common_factors, logfile)
    is_vv_identity_parallel(left, right, d, W=None, prefix_length, ncpus, verbose, logfile)

Where the arguments are as follows:

 - `left` and `right` is the identity to check

 - `d` is the dimension

 - `W` is an optional set of finite words

 - `prefix` is an optional prefix that turns `is_XX_identity` into a partial stest (this is indirectly used in the parallel versions)

 - `prefix_length` is used to chunk the subwords into different jobs: if it is
   set to `k` then there will be `2^k` jobs which correspond to the `2^k`
   possible prefixes.

 - the argument `logfile` can be set to either `sys.stdout` or a string in
   which case details about the execution will be written in a file with that
   name.

You can also list identities (as tuple over `{0,1}`) with:

- `def sv_identities(n, d, u_start, u_stop, nb_mats)`: where

  - `n` -- length of the identity

  - `d` -- dimension

  - `u_start` -- an optional start

  - `u_stop` -- an optional stop

  - `nb_mats` -- the number of matrices used to numerical rejection

Here are some examples

	sage: from max_plus import *

    sage: p = 'xyxxyy'
    sage: s = 'xxyyxy'
    sage: is_sv_identity(p+'x'+s, p+'y'+s, 3)
    True

    sage: p = 'xyyxxxyyyyxxxx'
    sage: s = 'yyyyxxxxyyyxxy'
    sage: %time is_sv_identity(p+'x'+s, p+'y'+s, 5)
	CPU times: user 48 ms, sys: 12 ms, total: 60 ms
	Wall time: 55.7 ms
	True

	sage: p = 'xyyxxxyyyyxxxxxyyyyy'
	sage: s = 'xxxxxyyyyyxxxxyyyxxy'
	sage: %time is_sv_identity(p+'x'+s, p+'y'+s, 6)
	CPU times: user 276 ms, sys: 8 ms, total: 284 ms
	Wall time: 259 ms
	True
	True
    sage: is_sv_identity(p+'x'+s, p+'y'+s, 7)
    False

	sage: p,s = vincent_sv_prefix_suffix(7)
	sage: is_sv_identity_parallel(p+'x'+s, p+'y'+s, 7, prefix_length=3, logfile=sys.stdout)
	u = xxxxyx
	num ext. occ.: 371
	num faces    : 16
	num verts    : 58
	polytope computation in 0.0998837947845secs

	u = xxxyxx
	num ext. occ.: 289
	num faces    : 11
	num verts    : 30
	polytope computation in 0.0397710800171secs
	...
	u = yxyyyy
	num ext. occ.: 371
	num faces    : 16
	num verts    : 58
	polytope computation in 0.0290420055389secs

	True
	sage: W = FiniteWords([0,1])
	sage: p = W((0,0,1,1,0))
    sage: s = W((0,0,1,0,1))
    sage: u = p + W([0]) + s
	sage: v = p + W([1]) + s
    sage: is_sv_identity_parallel(u, v, 3, W=W, prefix_length=2)
	True
    sage: t = WordMorphism({0:[0,1], 1:[1,0]}, domain=W)
    sage: is_vv_identity_parallel(t(u), t(v), 3, W=W, prefix_length=2)
	True

    sage: for i in sv_identities(11, 3):
    ....:     print i
	((0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0), (0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0))
	((0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0), (0, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0))
	((0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0), (0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0))
	((0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0), (0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0))
	((1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0), (1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0))
	((1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0), (1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0))
	((1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1), (1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1))
	((1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0), (1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0))
	((1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1), (1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1))
	((1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0), (1, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0))
	((1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0), (1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0))
	((1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0), (1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0))
	((1, 1, 0, 0, 1, 1, 0, 1, 1, 0, 0), (1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0))
	((1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1), (1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1))
	((1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1), (1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1))
	((1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0), (1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0))

## Contact

Any comment or remark is welcome at vincentDOTdelecroixATlabriDOTfr
