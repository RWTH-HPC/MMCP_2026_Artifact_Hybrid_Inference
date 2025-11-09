# Solver # {#dvSolver}

One of the central building blocks of our framework is the __solver__. Each solver represents one specific solution method, e.g. Discontinuous Galerkin, Finite Volume or Lattice Boltzmann, which can be used to solve a variety of partial differential equations.
The solver holds the solution state as well as auxiliary variables for each cell, which is typically organized in the corresponding [collector](@ref dvCollector).
Following functionalities __have__ to be provided by the solver:

* Initialization of data structures and solution state
* Solution algorithm
* Handle and apply boundary conditions
* Save solution state
* Save/load state to perform restart

A key concept for m-AIA is the joint grid for potentially multiple solvers. Thus, the solver itself does no hold any grid information, but a reference to the joint grid.
To be more specific, each solver owns a grid proxy which provides only the relevant information for the specific solver.
To learn more about the grid proxy concept, read [here](@ref dvProxy).

For some solvers you can chose between multiple system of PDEs. This is implemented via meta programming, i.e. these solvers have an additional template parameter which has to provide the equation specific parts.
This allows for compile-time optimization of the solution algorithm itself while maintaining the flexibility of choosing different systems of PDEs.

Multiple solvers can be coupled to solve a system of PDEs. This can be in the form of volume or surface coupling and is implemented via the [coupler concept](@ref dvCoupler).
Note, that solvers can't access each others variables directly, this is only done by the coupler, which can read/write from/to its participating solvers.

To support more advanced features like [adaptive mesh refinement](@ref dvAmr) and [dynamic load balancing](@ref lb), the solver needs to provide additional functionalities:

__AMR:__
* Define solution state based sensor values, used to determine whether a cell needs to be refined or coarsened
* Handle solution state for refined/coarsened cells

__DLB:__
* Define relevant cell types which should be taken into account during DLB
* Pack/unpack solution state during DLB
