# Level-set (LS) # {#ugLS}

[TOC]

The Level-set (LS) solver of m-AIA is used to track iso-surfaces. That is, the motion of a
rigid body, the deformation of a solid body, or contour surfaces of a flame front are tracked.
Although the LS solver can be used as an independent solver, it is coupled to other solvers in
most applications. In the following, the setup of two main applications are described. These
are rigid body motion and combustion simulations. Before these specific use cases are explained,
general properties are discussed.

## General adjustment and properties of the LS solver
To activate the level-set solver set
```
solvertype.* = "MAIA_LEVELSET_SOLVER"
```
When running a LS simulation several properties should be set in general. These are explained
in the tables below.

### G-Band Method
The level-set function \f$G\f$ is a signed distance function to a distinct surface. As described above, this can be a flame front or the surface of a rigid body. On the surface \f$G = 0\f$.
For a fast computation of the level-set, the \f$G\f$ field is only calculated in a band close to the level-set surface. Outside of this band a constant outside \f$G\f$ value is set. The width of the level-set band is set by
```
gBandWidth = 4
```
Here, 4 describes the number of cells on both sides of the interface.  
If the underlying Cartesian mesh is non-uniform, the level-set grid is only refined close the \f$G = 0\f$ contour, since no level-set value is computed outside the level-set band anyways.
The refinement is adjusted by the following two properties.  
The property `gShadowWidth` defines the distance from the \f$G = 0\f$ surface in which the level-set grid is refined. It is set in the same way as the `gBandWidth`.
Because the rigid bodies are moving, the refined grid region is continuously adjusted according to the body movement. Therefore, always set `adaptation = true`, if your grid is non-uniform.
The property `gInnerBound` describes at which point an adaptation of the grid is necessary. It describes a boundary inside the refined region and once a cell inside the level-set band intersects this boundary, the adaptation is forced.
The position of this boundary is `gShadowWidth-gInnerBound`.
For an efficient simulation ensure that `gShadowWidth-gInnerBound > gBandWidth + 2`.  

```
gShadowWidth = 10
gInnerBound = 2
```

## Rigid body motion
To activate the rigid body motion of the LS solver the `solverMethod` should be set to one of the
following options:  ??  

In the following, it is described how the LS solver is set up for a coupled LS FVMB simulation.
For this  
`solverMethod = "MAIA_RUNGE_KUTTA_MB_SEMI_LAGRANGE_LEVELSET"`
To allow the LS solver to run set
```
levelset = true
```
To inform the LS solver that it is coupled to a FVMB solver set
```
levelSetMb = true (in case the FVMB solver is coupled to the LS)
lb = true (in case the LB solver is coupled to the LS)
```  
The semi langrange level-set can either handle translatory or combined translatory and rotational motion. This is defined by the property `levelSetDiscretizationScheme`.
```
levelSetDiscretizationScheme = "BACKWARDS_PAR" # translatory motion
levelSetDiscretizationScheme = "ROTATING_LS" # translatory + rotational motion
```
In the following the combined translatory and rotational motion `levelSetDiscretizationScheme = "ROTATING_LS"` is described.
To use this discretization scheme set
```
LsRotate = true
```
To inform the level-set solver that it describes a forced body motion set
```
constructGField = false
trackMovingBndry = true
```
The level-set field can either be build from a specified analytical function or from a STL-file in which the surface of the \f$G = 0\f$ contour is defined. To activate the level-set reconstruction from an STL-file set
```
GFieldInitFromSTL = true # Build G-Field from STL
GFieldFromSTLInitMode = 1 # Calc G-Field from STL
initFromRestartFile = false # true does not work
```
The STL-files the level-set is reconstructed from need to be given in the geometry.toml file. The assignment of the STL-file to a specific level-set is handled by the property `bodyBndryCndIds`.
In the geometry.toml an id number is assigned to each STL-file
```
filename.0 = "stl/box.stl"
filename.1 = "stl/sphere.stl"
BC.0 = 3006
BC.0 = 3014
```
By entering these ids, the assignment is completed. Here, to STL files are assigned.
```
bodyBndryCndIds = [3006,3014] # STL from geometry.toml
```
The level-set solver can contain more than one level-set. The number of level-set functions is set by
```
maxNoLevelSets = 3
```
If more than one level-set function exists, a combined level-set function should be created to simply the coupling with the FVMB solver. Note that this combined level-set counts as an individual level-set. To create the combined level-set set
```
buildCollectedLevelSetFunction = true
```
To address the individual level-set functions in the code. A body to set table is generated. The properties below, inform the level-set solver that there are two level-set functions (Here the combined level-set is not counted). Each level-set is constructed from a single STL file. Here, two level-set functions are created from two STL files.
```
nodifferentSets = 2 # Combined level-set is not counted
bodiesinSet = [1,1] # Combined level-set is not counted
```
The algorithm reconstructing the level-set functions requires some user input. First, a point located inside the STL contour and inside the computational domain is required. Second, it must be defined wether the inside or the outside volume formed by the STL contour should have the positive level-set value. Those cells with a positive level-set value are considered for the flow simulation in the FVMB solver. In the example below two level-set functions exist.
```
initialInsidePoints = [-1.5, 6.7, 6.7, 0.0, 6.3, 6.3] # Combined level-set is not considered
levelSetSign = [1,-1,1] # Combined level-set is first value
```
In case a level-set is stationary the `computeSet` property can be set to false. This reduces the computational effort. The combined level-set function is considered by this property.
```
computeSet                     = [true,false,true] # Combined level-set is considered
```
The properties below define the movement of the rigid body. Several `bodyMovementFunctions` are implemented. Here, for both rigid bodies the case 8 is used. This is a translation of the body. The `amplitude` defines the velocity. Since we only want to consider a rotation of the body, `amplitudes` is set to zero for both bodies.
The last two properties define the rotation of the bodies. The rotation is always regarding the center of the coordinate system. `MaRot` sets the rotational speed and `bodyRadius` is a radius. Here, the second body rotates around the x-axis.
```
bodyMovementFunctions = [8,8]
amplitudes = [0.0, 0.0]
MaRot = [0.0, 0.0, 0.0, 0.0, 0.0,-0.1]
bodyRadius = [9.0, 10.0]
```
THe two properties below are required by the FVMB solver and need to be set in the given way.
```
### Required by fvmb
adaptLevelSetExtensionScheme = 2
LsMovement                   = true
```
The properties below need to be defined because no default value is set. However, since the reinitialization is not used. These properties are also unused.
```
### These properties need to be defined because no default value is set.
gRKMethod         = 5
extVelConvergence = 0.00001

### Reinitialization
reinitMethod      = "NONE"
gReinitIterations = 100
reinitConvergence = 0.0000000001
reinitCFL         = 0.05
guaranteeReinit   = false  # false for larger applications with higher refinement!
```

In the following it is described how the LS solver is set up for a coupled LS FVMB simulation.
For this  
`solverMethod = "MAIA_SEMI_LAGRANGE_LEVELSET_LB"`

`solverMethod = "MAIA_SEMI_LAGRANGE_LEVELSET"`

## Numerical time integration

As is the case of the FV solver, a Runge-Kutta method is used for the numerical time integration of the discretized level-set equation. 

An example of the properties is as follows, where the chosen RKMethod is set to 1, which activates a `three step third order TVD Runge-kutta scheme`. Therefore the `nogRKSteps`is also set to three to be consistent with the formulation of the method and the RK coefficients are chosen via the property `grkalpha-step`.

```
 gRKMethod = 1
 nogRKSteps = 3
 grkalpha-step = [1.0, 0.25, 0.6666666666666666]
```

## Solving the Jacobi-Hamilton equation

<table>
<tr><th>Property                         </th><th>Type  </th> <th>Description</th></tr>
<tr><td>levelSetDiscretizationScheme          </td><td> String   </td> <td> This property sets the level set discretization scheme of the level set equation.  \n
    possible values are:
    <ul>
    <li>US1 - first order upwind scheme</li>
    <li>UC3 - third order upwind scheme</li>
    <li>UC3_SB - third order upwind scheme</li>
    <li>UC5 - fifth order upwind scheme</li>
    <li>UC5_SB - fifth order upwind scheme</li>
    <li>WENO5 - fifth order upwind scheme</li>
    <li>WENO5_SB - fifth order upwind scheme</li>
    <li>UC11 - higher order upwind scheme</li>
    </ul></td></tr>
</table>


## Combustion simulations 

Should a combustion simulation with the level-set method be desired, the solverMethod is set to

`solverMethod = "MAIA_RUNGE_KUTTA_GEQU_PV"` .

### Reinitialization Method
Reinitialization of the level-set is necessary for combustion simulations using the level-set method. Solving the transport equation moves the zero level set \f$ \varphi_0 \f$ correctly, but may perturb the level-set function near \f$ \varphi_0 \f$, leading to large or small gradients of the level-set variable near \f$ \varphi_0 \f$. The property `reinitMethod` controls the type of reinitialization method used. A more thorough explanation of some of the methods can he found [here](@ref ).

All the convergence criteria and reinitialization control can also be chosen by the relevant properties. A concise summary of the most relevant ones can be found in the table below, however it is not only limited to these properties should a finer control be desired.

<table>
<tr><th>Property                         </th><th>Type  </th> <th>Description</th></tr>
<tr><td> reinitMethod          </td><td> String   </td> <td> This property sets the reinitialization method. \n
    Possible values are:
    <ul>
    <li>CR1 </li>
    <li>CR2 </li>
    <li>HCR1 </li>
    <li>HCR2 </li>
    <li>HCR2_LIMITED </li>
    <li>HCR2_FULLREINIT </li>
    <li>SUS5CR </li>
    <li>CR2PLUS </li>
    <li>RSU </li>
    <li>DL1 </li>
    <li>DL2 </li>
    <li>SUS_1 </li>
    <li>SUS_1PLUS </li>
    <li>SUS_2 </li>
    <li>SUS_WENO5 </li>
    <li>SUS_WENO5PLUS </li>
    </ul> </td></tr>
<tr><td> gReinitIterations          </td><td>Int   </td> <td> Triggers the maximum number of reinitialization steps to reach the convergence criterion \ref
    reinitConvergence</td></tr>
<tr><td> minReinitializationSteps          </td><td>Int   </td> <td> This property triggers the minimum number of reinitialization steps to reach the convergence criterion \ref
    reinitConvergence </td></tr>
<tr><td> reinitCFL          </td><td> Float   </td> <td> Sets the reinitialization cfl number which is used to calculate the pseudo time step for the
    reinitialization  \f$ \Delta t_{pseudo} = CFL_{reinit} * \Delta x \f$, with  \f$ \Delta x \f$ being the smallest G
    cell width. </td></tr>
<tr><td> reinitConvergence          </td><td> Float   </td> <td> Sets the convergence criterion of the reinitialization. The criterion should be of the order of
    \f$10^{-3} - 10^{-10}\f$. </td></tr>
</table>

## Adaptive mesh refinement
If you want to learn how to enable AMR for you simulation, read [here](@ref ugAMR). The list of all available sensors for the LS solver can be found @subpage sensorsLS "here".
