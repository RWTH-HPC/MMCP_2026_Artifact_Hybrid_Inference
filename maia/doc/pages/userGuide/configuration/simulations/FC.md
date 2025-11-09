# Finite Cell (FC) # {#ugFC}

[TOC]

In the following the Finite Cell (FC) solver is briefly discussed. That is, a short overview
about how to start the solver, the general settings, and an overview about boundary condition
settings is given. Further information about the mathematical modelling, i.e., the principle of
virtual work (PVW) and the numerical method of the FC is presented here ([PVW](@ref mmFCPVW) and
[FC](@ref nmFC)).

## Overview
The finite cell (FC) solver can be used to simulate structural mechanic problems. The solver
works in 2D and 3D. However, the 2D applications have not yet been tested. To start the solver,
the following properties have to be specified in the property file:
```
gridGenerator = false

flowSolver = true

solvertype.* = "MAIA_FINITE_CELL"
```

@note So far, only linear-elastic materials can be simulated. Furthermore, the simulation is static,
i.e., no dynamic forces can be applied.

## General settings
At first, simulation specific properties are set. These are, the dimension of the application, the
p-refinement, i.e., the number of nodes per cell, and the depth of the subcell integration per
boundary.
```
dim = 3

polyDeg = 1

subCellLayerDepth = [0, 0, 0, 0, 2]
```

These setting result in a 3D simulation using a [p-refinement](@ref nmFCpRef) with \f$ (p+2)(p+2)(p+2) = 27 \f$
nodes per cell, a [subcell integration](@ref nmFCsubCell) depth of \f$0\f$ at boundaries \f$0 - 3\f$, and a subcell
integration depth of \f$2\f$ at boundary \f$4\f$. Next, the [material properties](@ref mmFCMaterial), i.e., the Young's
modulus and the Poisson's ratio, are set. For example use:
```
EModule = 100000.0

PoissonRatio = 0.3
```

The method to solve the linear system of equations is the BiCGStab method. The number of maximum
iterations and the precision of the solution is set by:
```
epsBiCG = 1e-12

noIterations = 10000
```

To get extra debug output use:
```
testRun = true
```

To calculate the resulting stresses from the strains use
```
saveStress = true
```

## Boundary condition settings
In the FC solver, two different types of boundaries, i.e., natural and essential boundaries,
exist. At essential boundaries the displacement is set. At natural boundaries the traction
is set. The ids of the boundaries, to be specified in the geometry file of the simulation,
are given in [BC](@ref nmFCBC). Additionally, further properties are required in the properties
file of the simulation. These are the number of essential and natural boundaries, the
direction of the displacement and the load, and the load vector:
```
noFixedSeg = 2

noForceSeg = 2

segFixationDirs_0 = [1.0, 0.0, 0.0]

segFixationDirs_4 = [1.0, 1.0, 1.0]

segLoadDirs_1 = [1.0, 0.0, 0.0]

segLoadDirs_2 = [0.0, 1.0, 0.0]

segLoadValue_1 = [10.0, 0.0, 0.0]

segLoadValue_2 = [0.0, 10.0, 0.0]
```

In this setup, two essential and two natural boundaries exist. The displacement of the essential
boundaries is \f$ 0 \f$ in x-direction (boundary \f$ 0 \f$) and \f$ 0 \f$ in x-, y-, and
z-direction (boundary \f$ 4 \f$). At the natural boundaries, a traction of \f$ 10 \f$ is set
in x-direction (boundary \f$ 1 \f$) and y-direction (boundary \f$ 2 \f$).
