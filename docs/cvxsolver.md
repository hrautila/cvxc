
# NAME

cvxsolver - interior point conelp and convex solver

# SYNOPSIS

`cvxsolver -P filename [-v]`
`cvxsolver -c filename -G filename -h filename -d dimensions [-A filename -b filename -v]`
`cvxsolver -F name -c filename -G filename -h filename -d dimensions [-A filename -b filename -v]`
`cvxsolver -F name -G filename -h filename -d dimensions [-A filename -b filename -v]`

# DESCRIPTION

Runs convex solver on problem defined by input arguments. The solver is porting of Python CVXOPT solvers.

    minimize c*x^T
    s.t    G*x <= h
           A*x = b

See python package CVXopt documentation for more details.

# OPTIONS

-A --matrix-A
  Pathname of JSON serialized column major matrix.

-b --matrix-b
  Pathname of JSON serialized column major matrix.

-c --matrix-c
  Pathname of JSON serialized column major matrix.

-G --matrix-G
  Pathname of JSON serialized column major matrix.

-h --matrix-h
  Pathname of JSON serialized column major matrix.

-d --dims
  Pathname file or string holding JSON serialized problem dimensions. Dimensions data defines non-linear, 
  linear, second order cone and semidefinite constraint sizes. JSON serialization is an object with 
  following attributes:

  - nl: scalar, number of non-linear constrainst
  - l : scalar, number of linear constraints
  - nq: scalar, number second order cone constrainst
  - ns: scalar, number of semidefinite constraints
  -  q: array, SOCP constraint sizes
  -  s: array, SDP constraint sizes
      
  Currently the JSON serialization must have all scalara attributes preceding the array attributes. Attributes
  missing from serialization are assumed to be zero.

  -F --function
  Basename of shared object file that defines the convex constraint and convex target for CPL and CP problems.
  Shared object is loaded and initialization function *cp_init* is called with to arguments, first pointer
  to *cvxc_convex_program_t* structure and the second is opaque pointer to initialization data. The structure
  *cvxc_convex_program_t* is defined in 'cvxc.h' header.
  
  -P --program
  JSON serialization of complete problem parameter as defined in *cvxc_params_t*. See *cvxc.h* header
  
  -v
  Be somewhat verbose.

# CAVEATS

Missing documentation.

# AUTHOR

Harri Rautila, <harri.rautila@gmail.com>

# SEE ALSO

For details on solver parameter see Python CVXOPT documentation.



