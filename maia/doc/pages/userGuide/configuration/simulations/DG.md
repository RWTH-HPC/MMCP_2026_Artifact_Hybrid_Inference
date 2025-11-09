# Discontinuous Galerkin (DG) # {#ugDG}

[TOC]

The following guide is an exemplary recipe to setup a basic simulation with the DG solver.
For general settings (e.g. I/O), which are not specific to DG, see [general application control](@ref ugGeneralApplication). Properties in angle brackets
'''< >''' indicate all possible choices.
It is also encouraged to examine the testcases. The link to the repository can be found in the [references](#ugDG_ref).

# Problem Definition
Activate the DG solver and set the dimension of the spatial discretization:
```
  solvertype.default = "MAIA_UNIFIED"
  solvertype.0 = "MAIA_DISCONTINUOUS_GALERKIN"
  nDim = <2,3>
```
The DG solver supports space discretizations in 2D and 3D. The first step is to choose the right
system of equations to solve
```
  dgSystemEquations = <"DG_SYSEQN_ACOUSTICPERTURB", "DG_SYSEQN_LINEARSCALARADV">
```
For the properties specific to the available particular system of equations, continue
[here](@ref ugDG_sysEqn).
The space discretization by the DG method and the Runge-Kutta time integration can be controlled separately and are
outlined in the next section.

# Space Discretization # {#ugDG_SpaceDiscretization}
The solver operates on a Cartesian mesh and supports local grid refinements, which in the DG
terminology is called h-refinement. Additionally, the spatial resolution can be controlled by the
choice of the basis function, which for the collocation type approach also controls the numerical
quadrature
```
  dgIntegrationMethod = <"DG_INTEGRATE_GAUSS", "DG_INTEGRATE_GAUSS_LOBATTO">
```
As a rule of thumb Gauss DGSEM has superior dispersion properties (cf. [Gassner2011][Gassner2011]) and is preferred for acoustic
application. `DG_INTEGRATE_GAUSS_LOBATTO` satisfies the summation-by-parts property (SBP) and might
be superior in terms of numerical stability.
One advantage of the DG method is the ability to locally control the spatial discretization order
by adjusting the degree of the polynomial basis functions in each cell, called p-refinement
```
  initPolyDeg       = 3
  minPolyDeg        = 3
  maxPolyDeg        = 4
  pref              = 1                               // Activate p-refinement
  prefCoordinates_0 = [0.5, 0.0, 0.0, 1.0, 1.0, 0.75] // Box p-refinement [-x,-y,-z,+x,+y,+z]
  prefPolyDeg_0     = 4                               // Polynomial order of 0th patch
```
This example defines one p-refinement patch with index zero. Multiple p-refinement patches can be combined.
Note that the time step \f$\Delta t\f$ is a function of the polynomial degree of the basis functions, see [here](@ref nmDG).

# Time Discretization
The (temporal) integration method to be used in the DG solver can be set by
```
  dgTimeIntegrationScheme = "DG_TIMEINTEGRATION_CARPENTER_4_5" // 4th order, 5 stages
  cfl                     = 1.0
```
A complete list of all available integration schemes are listed here: @ref dgTimeIntegrationScheme.
The conversion from the `cfl` number to the time step \f$\Delta t\f$ is scaled for all combinations
of `dgSystemEquations` and basis functions to get for most applications the stability limit for `cfl = 1.0`.
Note, that in general it is recommended to adapt the spatial and temporal order of accuracy.
For coupled simulations an optimal interleaving of all solvers is achived with matching stage
numbers of the Runge Kutta schemes.


# System of Equations # {#ugDG_sysEqn}
The DG solver is designed for the solution of time dependent systems of first order convection
dominated PDEs. Besides the common configuration of the DG solver described in the previous
sections, depending on the underlying system of equations to be solved, further configuration is
required in particular for the initial and boundary conditions.

## APE
The acoustic perturbation equations in the [APE-4](@ref mmAPE) variant require the specification of
mean variables, so called `nodeVars`, which are time independent in addition to the initial acoustic
field, both set in DgSysEqnAcousticPerturb<nDim>::calcInitialCondition(). The type of initial
condition can be controlled by the property `initialCondition`. A list of all possible values can be
found @subpage dgICApe "here".

The acoustic waves can be excited by noise sources, which can be analytically described or taken
from a CFD solution. The property `sourceTerm` selects the source formulation. A list of all
possible values can be found @subpage dgSourceApe "here".

To fully define the physical problem proper boundary conditions must be specified in
`geometry.toml`. A list of all possible boundary conditions can be found @subpage dgBCApe "here".

### Standalone
The APE-4 system can be solved with prescribed noise sources for benchmarking and testing new DG features.

### Online Coupled
In the directly coupled setup the noise sources are exchanged in memory from a concurrently running
CFD simulation. The mean quantities and the source terms in the APE-4 are extracted from the respective CFD
simulation, and hence
```
  initialCondition  = 790 // nodeVars are not set & all fluctuations are set to zero
  sourceTerm        = 0
```
\todo Coupling related configuration

For more technical details about the coupling, see [here](@ref nmDgX).

## Linear advection equation
A list of all possible values for `initialCondition` can be found @subpage dgICLs "here".<br>
A list of all possible values for `sourceTerm` can be found @subpage dgSourceLs "here".

The theory to linear advection equations can be found [here](@ref mmLSA).


# References # {#ugDG_ref}
-  Details about the specific implementation of the DG scheme can be found [here](@ref nmDG) or in
[Schlottke2017][Schlottke2017].
-  A reference to the full list of properties can be found [here](#propertyPageDG).
-  A collection of testcases to start with can be found in the subversion repositories
[url](http://svn.aia.rwth-aachen.de/maia/testcases/DG/) for DG only and FV-DG coupled setups or
[url](http://svn.aia.rwth-aachen.de/maia/testcases/DG_LB/) for LB-DG coupled setups.

# Literature
-  Schlottke-Lakemper M. A Direct-Hybrid Method for Aeroacoustic Analysis. Dissertation, RWTH Aachen
University, 2017, [10.18154/RWTH-2017-04082][Schlottke2017].
[Schlottke2017]: http://publications.rwth-aachen.de/record/688887/files/688887.pdf
-  Gassner G, Kopriva D. A comparison of the dispersion and dissipation errors of Gauss and Gauss–Lobatto discontinuous Galerkin spectral element methods. SIAM Journal on Scientific Computing, 33(5), 2560-2579, 2011, [10.1137/100807211][Gassner2011].
[Gassner2011]: https://epubs.siam.org/doi/abs/10.1137/100807211 
