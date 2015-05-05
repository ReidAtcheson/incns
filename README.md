Overview
=========

Out of personal curiosity and interest I am writing a solver for the
incompressible Navier-Stokes equations on a fully
periodic domain. This repository contains my current progress.


Details
-------
I start with a standard tensor-product Fourier basis
for each component of the fluid velocity, and I furthermore
combine this with a Galerkin approximation approach. Ideally
this should yield an exponential convergence rate and therefore
minimize memory requirements to simulate fine scale behavior;
however, handling pressure is tricky. Therefore to make
this a little easier I simply take the curl of the tensor-product
basis, this forces the basis to be divergence-free and eliminates
the pressure term. The curl-tensor-product basis yields a set
of functions which is highly linearly dependent, and so 
my initial fix to this was simply to use a rank-revealing QR
factorization supplied by Julia's "qrfact" function and then
extract the permutation vector. This yields both the numerical rank
of the basis as well as a reasonable set of basis functions to throw out
as redundant. That approach however severely limits the order I can use
because qrfact must be used on a dense matrix, and the dimensions
grow rapidly, quickly swamping out memory. What I am currently
attempting to do now is to detect the rank of the matrix while 
using the fact that the curl-tensor-product basis yields a very 
sparse mass matrix, allowing me to scale the approach to much
higher orders.

So far the Galerkin inner products are partially computed analytically
and partially computed through the Trapezoid rule, which is a rapidly
convergent quadrature rule for smooth periodic functions. Early
simulations simply construct the dense interpolation matrices
explicitly so as to make use of all cores through OpenBLAS. Later
work can use more efficient analytical formulas and some
threading model to utilize all cores.




