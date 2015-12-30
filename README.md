# A small Python library to deal with max-plus matrices relations

To run the library you need to install Sage (http://sagemath.org). Then
clone this repository (or download all the files). Go to the directory
where you downloaded the file. Start Sage and run

    sage: %runfile int_max_plus.pyx
    sage: %runfile max_plus.py

To create symbolic triangular matrices use the functions
`symbolic_max_plus_matrices_band(d, n, diag='v', surdiag='v', ch=None)` where

- `d` is the dimension

- `n` is the number of matrices

- `diag` describes whether the diagonal is constant equal to 0 (`'c'`), similar
  in each of the `n` matrices (`'s'`) or variable (`'v'`)

- `surdiag` describes whether the surdiagonal is constant, similar or variable.

- `ch`: the convex hull engine

For example you can do

    sage: x,y = symbolic_max_plus_matrices_band(3, 2, 'c', 'v')
    sage: x
    A 3x3 symbolic max plus matrix on 4 variables
    sage: print x.str()
	[  0 x_0 -oo]
	[-oo   0 x_1]
	[-oo -oo   0]
	sage: y
    A 3x3 symbolic max plus matrix on 4 variables
	sage: print y.str()
	[  0 x_2 -oo]
	[-oo   0 x_3]
	[-oo -oo   0]
    sage: p = x*y*x
    A 3x3 symbolic max plus matrix on 4 variables
    sage: print p.str()
	[   0 max(x_2, x_0) max(x_1+x_2, x_0+x_3, x_0+x_1) ]
	[ -oo             0                  max(x_3, x_1) ]
	[ -oo           -oo                              0 ]

From a symbolic matrix you can obtain an integer max plus matrix using the method `eval`

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

For full matrices you need to use the function `symbolic_max_plus_matrices(d,
n, ch=None, sym=False)` where

- `d` is the dimension

- `n` is the number of matrices

- `ch` is the convex hull engine

- `sym` specifies the implementation: if `sym` is `True` uses matrices that
  store only two coefficients (and deduce the others by symmetry). If `sym` is
  `False` uses matrices that store all coefficients.
