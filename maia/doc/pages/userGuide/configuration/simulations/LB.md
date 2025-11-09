# Lattice Boltzmann (LB) # {#ugLB}
[TOC]
This page gives an overview over the most important settings for the LB solver.  
If you are new to m-AIA make sure to first read the @ref ugGettingStarted section.  
For the general configuration of m-AIA please refer to @ref ugGeneralApplication.

# Example property file
<details>
<summary>Minimum properties_run.toml file from testcase LB/3D_sphere_Re100_parallel</summary>
```{.py}
# -------------- APPLICATION SETTINGS  --------------
 flowSolver                 = true
 restartFile                = false
 scratchSize                = 20.0
 outputDir                  = "out/"
 geometryInputFileName      = "geometry.toml"
 gridInputFileName          = "grid.Netcdf"
 timeSteps                  = 100
 solutionInterval           = 50
 restartInterval            = 5000  # must be an integer multiple of solutionInterval
 residualInterval           = 25
 noDomains                  = 4
# -------------- NUMERICAL PROPERTIES --------------
 executionRecipe            = "RECIPE_BASE"
 solvertype                 = "MAIA_LATTICE_BOLTZMANN"
 solverMethod               = "MAIA_LATTICE_BGK"
 nDim                       = 3
 noDistributions            = 19
 noSolvers                  = 1
# ------------ FLOW INITIALIZATION -----------
 initMethod                 = "LB_LAMINAR_INIT_PX"
# -------------- FLOW VARIABLES --------------
 Ma                         = 0.1
 Re                         = 100
 referenceLengthLB          = 16.0  # characteristic diameter in cell units on highest level
# ----------- GEOMETRY SEGMENTATION ----------
 inflowSegmentIds           = [0, 1, 2, 3, 4, 5] # id's of open-boundary segments
 initialVelocityVectors     = [1.0, 0.0, 0.0,
 			                   1.0, 0.0, 0.0,
			                   1.0, 0.0, 0.0,
			                   1.0, 0.0, 0.0,
			                   1.0, 0.0, 0.0,
			                   1.0, 0.0, 0.0] # velocity vectors for open boundaries
# -------------- GRID AND REFINEMENT PROPERTIES --------------
 maxNoCells                 = 150000  # collector size for flow solver
 reductionFactor            = 1.28  # the initial cell is length is factorized by this
 minLevel                   = 4
 maxRfnmntLvl               = 9
 interpolationDistMethod    = "sphere"  # method for measuring of wall distances
```
</details>

# Example geometry file
<details>
<summary>Minimum geometry.toml file from testcase LB/3D_sphere_Re100_parallel</summary>
```{.py}
 noSegments = 7
 inOutSegmentsIds = [0,1,2,3,4,5]

 filename.0 = "stl/back.stl"
 filename.1 = "stl/front.stl"
 filename.2 = "stl/top.stl"
 filename.3 = "stl/bottom.stl"
 filename.4 = "stl/inflow.stl"
 filename.5 = "stl/outflow.stl"
 filename.6 = "stl/sphere.stl"

 "body_segments.side1" = 0
 "body_segments.side2" = 1
 "body_segments.side3" = 2
 "body_segments.side4" = 3
 "body_segments.in" = 4
 "body_segments.out" = 5
 "body_segments.body" = 6

 BC.0 = 4030
 BC.1 = 4030
 BC.2 = 4030
 BC.3 = 4030
 BC.4 = 1000
 BC.5 = 4000
 BC.6 = 2000
```
</details>

# Numerical methods
### solverMethod
@par Type
  `string`
@par Default
  `none, required`

The property `solverMethod` determines which collision step is used in the LB solver.  
A List of the available collision steps and their respective values of the property `solverMethod` is generated here:
@subpage collisionStep.  
@note Some collision steps can only be used with a certain number of dimensions or distributions.  

<!-- ### interpolationType -->
### interpolationDistMethod
@par Type
  `string`
@par Default
  `perpOp`

The property `interpolationDistMethod` controls how the distance of cells to walls (defined by STLs) is calculated.
Possible values are:
- `STD` Should probably the standard method for 3D. Most accurate, but slow for STLs with a large number of triangles.
Not available for 2D.
- `perpOp` Uses the perpendicular operator.
- `sphere` Analytical solution for a sphere of radius 0.5 around the origin.
- `pipe` Analytical solution for a cylinder through the origin.

The different versions are implemented in the function LbBndCndDxQy::calculateWallDistances().

@note At the moment `perpOp` remains as the default value, but it is recommended to the more accurate `STD` option for
3D simulations.

### interfaceMethod
@par Type
  `string`
@par Default
  `FILIPPOVA`

Selects interface method for locally refined meshes.  
- `FILIPPOVA` Spatial and temporal interpolation of variables, rescaling of non-eq components
- `ROHDE` Volumetric formulation (no interpolation, no rescaling)
# General settings
### Re
@par Type
  `float`
@par Default
  `none, required`

The Reynolds number for the LB solver.
### Ma
@par Type
  `float`
@par Default
  `none, required`

The Mach number for the LB solver.
## Reference length
There are three different (but exclusive) properties that can be used to specify the reference length for the LB solver.
The reference length is the length in the Reynolds number.
@note Exactly one of these properties has to be present for the LB solver to run.
### referenceLength
@par Type
  `float`
@par Default
  `none`

This property sets the reference length in stl units, whereas internally the reference length is converted to LB units.
### referenceLengthLB
@par Type
  `float`
@par Default
  `none`

This property sets the reference length in multiples of the cell length at the finest refinement level. This corresponds
to the reference length that is internally used for nondimensionalization in the LB solver.
### referenceLengthSegId
@par Type
  `int`
@par Default
  `none`

With this property a segment Id (as used in the geometry file) can be specified that is used for the calculation of the reference lenght. For this
segment the hydraulic diameter is calculated. The method is described in more detail in the documentation for the
inline function LbSolver::calcCharLenAll(MInt segmentId).

## Forces
## Other settings
### updateAfterPropagation
@par Type
  `boolean`
@par Default
  `false`

Determines if the variables are updated at the end of the time step (after propagation) or at the beginning of the next
time step (before the collision). For the default value `false` the variables that are written to solution files are
actually from the previous time step, since they have not been updated. This can also cause problems for coupling
purposes
@note For now only implemented for the @ref sec_coll_cumulant

# Initial Conditions
The LB solver provides several initialization conditions. Some of these are quite specific. Others are more general.
In the following, only the general ICs are presented. A full list of ICs can be found here: @subpage LBICList

## initMethod
@par Type
  `string`
@par Default
  `none, required`

This property determines the initialization procedure during which the initial state (variables and distributions) is
set for all cells. The following list documents only the most common methods. For a complete list refer to the
documentation of @ref LBinitMethod.
- `LB_FROM_ZERO_INIT` All velocities set to zero, all other variables set to their infinity state. Distributions in 
equilibrium state.
- `LB_LAMINAR_INIT_PX` All `initMethod`s starting with `LB_LAMINAR_INIT_` set the velocity in one Cartesian direction
to the infinity state, while all other velocities remain zero. The last charachter `X`, `Y`, or `Z` determines the 
Carteisan direction. The preceding `P` (plus) or `M` (minus) determines the sign of the velocity.
- `LB_LAMINAR_INIT_MX` see above
- `LB_LAMINAR_INIT_PY` see above
- `LB_LAMINAR_INIT_MY` see above
- `LB_LAMINAR_INIT_PZ` see above
- `LB_LAMINAR_INIT_MZ` see above

@note Some initialization methods that generate synthetic turbulent flow require additional properties,
e.g. @ref FFTInit and @ref noPeakModes

## Ramped initialization
To prevent large shocks or instabilities due to initialization with a high Reynolds number, ramping of the Reynolds
number periodicCartesianDir = 1, 0, 0is supported. If ramped initialization is enabled by the @ref tanhInit property, the simulation starts at a lower
Reynolds number defined by @ref initRe. After a number of timesteps that is defined by @ref initStartTime the Reynolds
number is increased following a tanh function. The final Reynolds number (defined by the property `Re`) is reached after
additional @ref initTime time steps.

# Boundary Conditions
The boundary conditions have to be specified for all boundary segments in the `geometry.toml` file.  
A complete list of all boundary conditions is generated here: @subpage LBBCList

There are no strict rules about the naming scheme of the boundary conditions. But as a rule of thumb you can expect the
following types of boundary conditions for these Ids:
- `1xxx` Velocity BCs (e.g. inflow)
- `2xxx` Bounce back BCs (e.g. walls)
- `4xxx` Pressure BCs (e.g. outflow)

## Segmend Ids
A segment ID is assigned to each boundary segment in the `geometry.toml` file. 

### inOutSegmentIds
In the `geometry.toml` file it has to be defined,
what segment ID is assigned to either an inflow or outflow boundary.
For such a segment, an extrapolation direction is required,
since some of the inflow and outflow boundary conditions require extrapolation from the inner cells
to calculate flow quantities at the boundary cells.

The extrapolation directions are defined with the property `initVelocityMethod` in the `properties_run.toml` file.
If the property is set to `read`, the normal direction of the boundary segment needs to be specified with an array
defined by the property `initialVelocityVectors`.
With the option `calcNormal` the normal direction is calculated from vertices of the boundary segment,
The option `fromSTL` reads the normal direction from an STL file.

### inflowSegmentIds
With the property `lbControlInflow` in the `properties_run.toml` file,
the velocity profile at inflow boundaries can be prescribed.
Per default (`lbControlInflow=0`) a profile with a uniform inflow velocity is used. 
If an inflow boundary shall have a parabolic velocity profile (`lbControlInflow=1`),
its segment ID needs to be added to the property `inflowSegmentIds` in the `properties_run.toml` file.

### periodicSegmentsIds
For imposing periodic boundary conditions the corresponding segment IDs need to be specified
in the `geometry.toml` file with the property `periodicSegmentsIds`.
Additionally, the property `periodicCartesianDir` needs to be determined in the `properties_run.toml` file,
in form of a vector that points into the direction in which the grid should be periodic.


## Moving Boundary

## Adaptive mesh refinement
If you want to learn how to enable AMR for you simulation, read [here](@ref ugAMR). The list of all available sensors for the LB solver can be found @subpage sensorsLB "here".
