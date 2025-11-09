# Finite Volume Moving-Boundary (FVMB) # {#ugFVMB}

## Overview
The finite-volume moving-boundary (FVMB) solver enhances the regular finite-volume ([FV](@ref ugFV)) solver by an immersed-boundary method. In this, the FVMB solver inherits from the FV solver and is thus able to access all of its parent variables and functions. The general functionality is shared ([ER](@ref exR)) whereas specific functions are adapted to accomodate for the moving bodies.

The implemented immersed-boundary method allows for the embedment of rigid bodies of various shapes and properties. The FVMB solver can therefore be used to model a wide variety of environments such as rigid particles in a particle-laden flow or complex geometrical settings as in, e.g., combustion engines with moving pistons.
The information used to model the geometric shape is provided by the Levelset solver ([LS](@ref ugLS)) whereas particles of varying shape can be directly modeled by an analytical levelset (ANAL-LS) incorporated into the FVMB solver. In this, the term *levelset* refers to a distance-based function \f$ \phi \f$ where \f$ \phi>0 \f$ comprises to the outside, \f$ \phi<0 \f$ to the inside of a body, and \f$ \phi=0 \f$ to the body surface. 

### Fluid-Body interaction
To discretize the body based on the given levelset information, a Cartesian multi cut-cell method is used ([CutCell](@ref nmFVCC)). All conservative units, i.e., mass, momentum, and energy, are strictly conserved at the fluid-particle interfaces 

### Surfaces
Separate from the source of information (LS vs ANAL-LS), all boundaries of the body are represented by boundary surfaces. In this, the aforementioned multi-cut cell algorithm triangulates the levelset information and adds surfaces where appropriate. Each surface is then supplied a boundary ghost cell to describe the properties at the body surface ([GCB](@ref nmFVCC)). Multi-cut cells, that is, cells with multiple boundary surfaces, are assigned multiple ghost cells.

### Boundary conditions
Different boundary conditions can be enforced onto the moving-boundary surfaces. Commonly used options include adiabatic boundary conditions (3006), isothermal boundary conditions (3008), and boundary conditions for static (non-moving) boundaries (3600). All boundary conditions are specific for (non-)moving boundaries and can be applied to the Navier-Stokes- as well as the Euler-equations.



### Movement
The immersed boundaries are updated according to different properties/algorithms and can be classified into two different classes: forced movement and free motion. 
Forced movements are commonly used by complex geometrical boundaries such as multi-stage axial turbines and are computed in the LS solver and later on transferred into the FVMB solver ([Coupling](@ref nmCoupling)). Here, the FVMB solver receives updated levelset information used for the discretization of the moving boundary surfaces whereas the movement update is compute in different solvers.
The concept of free motion is generally described by solving Newtons law of motion, time intregration of the resulting forces, and is commonly applied in particle computations in the FVMB solver directly. Here, each particle motion, that is, e.g., rotational and translational acceleration and velocity, is described by the integrated forces acting on its boundary reduced to its center of gravity.


The utilized methods and the holistic approach allows the direct investigation of the direct transfer of kinetic energy via the flux of linear and angular momentum between the fluid and the rigid body.
The fluid-boundary computations are further stabilized by a correction of intersected fluid cells whichs volume falls below a specified threshold.

The overall method exhibits very good stability and accuracy properties even for large displacements of the boundaries. 


## Runtime
A five-step Runge-Kutta Scheme of second-order accuracy is commonly applied to the temporal integration of the conservative variables while using the Euler equations as governing equations is possible. The computation of the RK Scheme is equivalent to that of the FV Solver ([FV](@ref ugFV)). As mentioned in the introduction, the FVMB Solver inherits from the FV Solver and thus, the general computing procedure is similar and, in some places, equal to its parent solver. A general overview will be given here, major differences are highlighted and similarities between both solvers will be omitted. For further reference, the reader is reffered to ([FV](@ref ugFV)).

#### Execution Recipe
The execution recipe describes the general workflow of the solver, different execution recipes are possible ([ER](@ref exR)). The main objective here is that different modalities for the update of the moving boundary bodies are possible. For example, the modelling of rigid particles in particle-laden flow is approached using a Predictor-Corrector Scheme which runs multiple cycles within a singular timestep which is in contrast to the investigation of geometries using forced motion such as a cylinder in a combustion engine.

All execution recipes have in common, that the call solversteps in a predefined order. A general overview of the used functions by the FVMB solver are given in the following.

#### PreTimestep
Marks a new timestep and, depending on the simulation settings, advances the moving boundary bodies properties and computes the resulting levelset information.

#### PreSolutionstep
Determines cell settings and marks moving boundary cells, generates the specific moving boundary cells based on the cell levelset values, computes the relevant boundary surfaces and adds ghost cells accordingly. The solver then updates the cell properties due to the moving of the surfaces based on, e.g., the redistribution of fluid mass inside the domain and in close proximity to the body boundaries.

#### Timestep
Integrates the governing equations in time and subsequently performs the relevant actions on the left-hand side and the right-hand side of the equations.

#### PostSolutionstep
Performs computations after the timestep integration. Herein, the residual is checked and, depending on the executionrecipe, the correction of the rigid body motion parameters and levelset information is updated based on the resulting boundary surface forces.

#### PostTimestep
Finalizes the timestep.


## Properties
The FVMB Solver has the majority of properties in common with the FV solver and thus, only major moving-boundary specific properties will be listed and described shortly in the following.

### Timestep
<table>
<tr><th>Property                         </th><th>Type  </th> <th>Description</th></tr>
<tr><td>timeStepAdaptationStart          </td><td>Int   </td> <td>Start of adaptation of the timestep value after this timestep</td></tr>
<tr><td>timeStepAdaptationEnd            </td><td>Int   </td> <td>The timestep at which a desired cfl number cflTarget is reached</td></tr>
<tr><td>cflInitial                       </td><td>Float </td> <td>The initial CFL number for the time step adaptation</td></tr>
</table>


### Moving-Boundaries
<table>
<tr><th>Property                         </th><th>Type    </th> <th>Description                                 </th></tr>
<tr><td>bodySamplingInterval             </td><td>Boolean </td> <td>Defines the sampling interval of the immersed bodies</td></tr>
<tr><td>trackBodySurfaceData             </td><td>Boolean </td> <td>Triggers the computation of relevant variables at the moving boundary surfaces and determines the resulting body forces </td></tr>
<tr><td>trackMovingBndry                 </td><td>Boolean </td> <td>Triggers the displacement of bodies in the moving boundary solver using the levelset field </td></tr>
<tr><td>trackMbStart                     </td><td>Int     </td> <td>For time steps smaller than the specified value, the bodies are not displaced and the levelset field is not updated</td></tr>
<tr><td>trackMbEnd                       </td><td>Int     </td> <td>For time steps larger than the specified value, the bodies are not displaced and the levelset field is not updated </td></tr>
<tr><td>complexBoundaryForMb             </td><td>Boolean </td> <td>Triggers if complex boundaries should be considered in the FVMB solver, i.e., cut cells are based on information from multiple levelset functions </td></tr>
<tr><td>generateOuterBndryCells          </td><td>Boolean </td> <td>Controls whether outer (non-moving) boundary cells should be created </td></tr>
<tr><td>gravity                          </td><td>Float   </td> <td>Gravitational acceleration in X/Y/Z direction </td></tr>
<tr><td>motionEquation                   </td><td>Int     </td> <td> Specifies the equations that are used for updating the body properties such as velocity. Possible values are:  
0: forced motion  
1: free motion under gravity  
2: Ahn & Kallinderis scheme  
3: Borazjani scheme  
4: Adaptation by Schneiders  </td></tr>
<tr><td>movingBndryCndId                 </td><td>Int   </td> <td>Specifies the boundary condition id assigned to the moving boundary surfaces </td></tr>
<tr><td>noEmbeddedBodies                 </td><td>Int   </td> <td>Number of distinct bodies </td></tr>
<tr><td>bodyTypeMb                       </td><td>Int   </td> <td>Select the moving boundary body type to simulate. Possible values are:  
0: generic body (STL)  
1: sphere  
2: piston  
3: ellipsoid  
4: NACA00XX  
5: split sphere  
6: cube  
7: tetrahedron  </td></tr>
<tr><td>bodyRadii                         </td><td>Float  </td> <td>Specifies the bodyradii used for the generation of spherical and ellipsoidal particles </td></tr>
<tr><td>densityRatio                      </td><td>Float  </td> <td>Sets the density of a moving body based on rhoInfinty </td></tr>
<tr><td>initialBodyCenter                 </td><td>Float  </td> <td> Specifies the initial body position of the particles</td></tr>
<tr><td>initialBodyVelocity               </td><td>Float  </td> <td>Specifies the initial body velocity of the particles </td></tr>
<tr><td>initialBodyRotation               </td><td>Float  </td> <td>Specifies the initial body rotation of the particles</td></tr>
<tr><td>fixedBodyComponents               </td><td>Int    </td> <td>Restricts translational movement of the rigid bodies </td></tr>
<tr><td>fixedBodyComponentsRotation       </td><td>Int    </td> <td>Restricts rotational movement of the rigid bodies </td></tr>
<tr><td>applyCollisionModel               </td><td>Boolean</td> <td>Enables the collision model (based on Glowinski et al.) </td></tr>
<tr><td>volumeLimitWall                   </td><td>Float  </td> <td>Fractional cell volume which triggers the Smallcell correction</td></tr>
<tr><td>volumeLimitOther                  </td><td>Float  </td> <td>Fractional cell volume which triggers the Smallcell correction</td></tr>
</table>

### Levelset
<table>
<tr><th>Property                          </th><th>Type    </th> <th>Description                                 </th></tr>
<tr><td>levelSetMb                        </td><td>Boolean </td> <td>Enables the use of levelset information for the computation</td></tr>
<tr><td>buildCollectedLevelSetFunction    </td><td>Boolean </td> <td>Triggers if the levelset functions from all sets are collected in a combined set </td></tr>
<tr><td>maxNoLevelSets                    </td><td>Int     </td> <td>Maximum number of levelsets used for the computation</td></tr>
<tr><td>constructGField                   </td><td>Boolean </td> <td>Enables the usage of levelset information for rigid particles </td></tr>
</table>

## Adaptive mesh refinement
If you want to learn how to enable AMR for you simulation, read [here](@ref ugAMR). The list of all available sensors for the FVMB solver can be found @subpage sensorsFVMB "here". Additionally, the sensors from the classic [FV solver](@ref sensorsFV) can be used.

## Example
In the following, a minimal example of a free falling, i.e., gravity driven, ellipsoidal particle in a rectangular periodic box is given. 
A single particle (*noEmbeddedBodies*) of ellipsoidal shape (*bodyTypeMb*,*bodyRadii*,*bodyRotation*) has an initial velocity of 0.1 in +x direction (*initialBodyVelocity*) and gets accelerated by gravitational forces (*gravity*). The bounding box uses cutoff boundary conditions which are additionally, in x-direction, periodic simulating a freefalling particle in an infinite length domain.

Beside the properties mentioned to compute the particle data, additional properties for different HPC relatived mechanisms are given , please refer to the respective pages of this documentation for further reading (ref general application,  fv, adaptation, balance, hpc, cutoff).

```
#============ STARTUP SETTINGS ============
gridGenerator                = false
flowSolver                   = true
noDomains                    = 24
restartFile                  = false


#============ DATA SETTINGS ============
noSolvers                     = 1
outputDir                     = "out/"
gridInputFileName             = "grid.Netcdf"
solutionOutput                = "./out/"
geometryInputFileName         = "geometry.toml"


#============ OUTPUT INTERVALS ============
solutionInterval              = 100
residualInterval              = 1
dragOutputInterval            = 1
restartInterval               = 250


#============ ADAPTATION ============
adaptation                     = true
allowInterfaceRefinement       = true
sensorAdaptation               = true

sensorType                     = ["INTERFACE","VORTICITY"]
sensorWeight                   = [-1.0,1.0]

adaptationInterval             = 50
adaptationStart                = -1


#============ BALANCE ============
balance                        = true
balanceAfterAdaptation         = true

loadBalancingMode              = 0
loadBalancingInterval          = -1


#============ SOLVER TYPE ============
solvertype.default             = "MAIA_UNIFIED"
solverMethod                   = "MAIA_RUNGE_KUTTA_MB_LEVELSET"
solvertype.0                   = "MAIA_FV_MB"
executionRecipe                = "RECIPE_ITERATION"
maxIterations                  = 2

solverOrder_0                  = 1
adaptationOrder                = 1

viscousFluxScheme              = "FIVE_POINT"
surfaceValueReconstruction     = "HOCD"
govEqs                         = "NAVIER_STOKES"


#============ TESTCASE SETTINGS ============
nDim                            = 3
initialCondition                = 0
Re                              = 100.0
Ma                              = 0.1
gamma                           = 1.4


#============ NUMERICAL SETTINGS ============
timeSteps                       = 3000

cfl                             = 1.0
cflViscous                      = 0.25

timeStepMethod                  = 100
upwindCoefficient               = 0.5
globalUpwindCoefficient         = 0.01

volumeLimitWall                 = 0.5
volumeLimitOther                = 0.5

noRKSteps                       = 5
rkalpha-step                    = [0.25, 0.1666666666666667, 0.375, 0.5, 1.0]


#============ GRID GENERATION ============
minLevel                        = 7
maxUniformRefinementLevel       = 7

maxBoundaryRfnLvl               = 8
maxRfnmntLvl                    = 8


#============ COLLECTOR SIZES ============
maxNoCells                       = 500000
maxNoSurfaces                    = 1000000
maxNoBndryCells                  = 50001

maxNoBndryCndIds                 = 3
scratchSize                      = 20.0


#============ MOVING BOUNDARY SETTINGS ============
trackMovingBndry                = true
trackBodySurfaceData            = true
complexBoundaryForMb            = true
generateOuterBndryCells         = true

constructGField                 = true
levelSetMb                      = true

maxNoLevelSets                  = 2

bodyFaceJoinMode                = 1
multiCutCell                    = 1

motionEquation                  = 1
bodyToSetMode                   = 10
buildCollectedLevelSetFunction  = true

bodyTypeMb                      = 3
bodyRadius                      = 2.0
bodyRadii                       = [2.0, 1.0, 1.0]
noEmbeddedBodies                = 1
densityRatio                    = 100.0

initialBodyCenter               = [ 0.0, 0.0, 0.0]
initialBodyRotation             = [ 25.0, 45.0, 65.0]
initialBodyVelocity             = [ 0.1, 0.0, 0.0]

fixedBodyComponents             = [0, 0, 0]
fixedBodyComponentsRotation     = [0, 0, 0]

movingBndryCndId                = 3006

adaptiveGravity                 = false
gravity                         = [1,0,0]

periodicDir                     = [1,0,0]


#============ CUTOFF SETTINGS ============
createBoundaryAtCutoff           = true

cutOff                           = 4
cutOffMethod                     = "P-P"

cutOffDirections                 = [0,1]
cutOffCoordinates                = [-40.0, 0.0, 0.0, 1.0, 0.0, 0.0, 40.0, 0.0, 0.0, -1.0, 0.0, 0.0]
cutOffBndryIds                   = [ 0, 0]
```
