# Sage library for max-plus matrix identities

To run the library you need to install Sage (http://sagemath.org). Then
clone this repository (or download all the files). Go to the directory
where you downloaded the file. Start Sage and run

    sage: %runfile int_max_plus.pyx
    sage: %runfile max_plus.py

Once that done, you can create band and upper triangular matrices
with the functions

    symbolic_max_plus_matrices_band(d, n, diag, surdiag)
	symbolic_max_plus_matrices_upper(d, n, diag, surdiag)

where

- `d` is the dimension

- `n` is the number of matrices

- `diag` describes whether the diagonal is constant equal to 0 (`'c'`), similar
  in each of the `n` matrices (`'s'`) or variable (`'v'`). This is optional and
  default to `'v'`.

- `surdiag` describes whether the surdiagonal is constant, similar or variable.
  This is also optional.

For example you can do

    sage: x,y = symbolic_max_plus_matrices_band(3, 2, 'c', 'v')
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

	sage: x2,y2,z2 = symbolic_max_plus_matrices_upper(3, 3)
    sage: x2.num_variables()
    18
	sage: print x2
	[  x0   x3  x4 ]
	[ -oo   x1  x5 ]
	[ -oo  -oo  x2 ]
	sage: print y2
	[  x6   x9  x10 ]
	[ -oo   x7  x11 ]
	[ -oo  -oo   x8 ]
	sage: print z2
	[ x12  x15  x16 ]
	[ -oo  x13  x17 ]
	[ -oo  -oo  x14 ]

From a symbolic matrix you can obtain an integer max plus matrix using the
method `eval`

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
    sage: pv = p.eval((1,-1,0,3))
    sage: pv
	[   0   1 4 ]
	[ -oo   0 3 ]
	[ -oo -oo 0 ]
    sage: pv == xv*yv*xv
    True

For full matrices you need to use the function

    symbolic_max_plus_matrices(d, n, ch=None, sym=False)

where

- `d` is the dimension

- `n` is the number of matrices

- `ch` is the convex hull engine

- `sym` specifies the implementation: if `sym` is `True` uses matrices that
  store only two coefficients (and deduce the others by symmetry). If `sym` is
  `False` uses matrices that store all coefficients.

For example

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


Any comment or remark is welcome at vincentDOTdelecroixATlabriDOTfr
