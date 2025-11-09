# Finite volume (FV) # {#ugFV}

## Overview
First, the flow solver has to be turned on, and the grid solver off by selecting following boolean values in the properties file: 

`gridGenerator = false`  

`flowSolver = true`

When starting an FV simulation, the user has to choose the following solver property :  

`solvertype.0 = "MAIA_FINITE_VOLUME"`

Next, an appropriate solver method for the time integration has to be chosen. For the finite volume solver the Runge-Kutta method can be chosen as:  

`solverMethod = "MAIA_RUNGE_KUTTA"`

The dimensionality of the simulation can be chosen by setting the value of `nDim` to 2 or 3. 

## Choosing the set of equations

The set of equations to be solved is implemented in the FV-solver through the sysEqn class. It can be set in the properties.toml file through the 
`fvSystemEquations` property.  

Following options exist:
* __Navier-Stokes__ (FV_SYSEQN_NS): Standard Navier-Stokes equations, non-dimensional, closed with the ideal gas equation. Constant heat capacities are assumed. 
* __RANS__ (FV_SYSEQN_RANS): Reynolds-Averaged-Navier-Stokes equations. Similar to the aforementioned NS-equations, but require of additional modelling and equations to close the system. 
* __Detailed Chemistry__ (FV_SYSEQN_DETCHEM): Dimensional Navier-Stokes equations expanded with additional transport and reaction equations. Non-constant heat capacities. 
* __EE-Gas__ (FV_SYSEQN_EEGAS): ...

<table>
<tr><th>Property                         </th><th>Type  </th> <th>Description</th></tr>
<tr><td> gamma         </td><td> Float </td> <td> Ratio of specific heats. Default value of 1.4 for fluid simulations concerning air.</td></tr>
</table>

#### Infinity values
Infinity values for the primitive variables are defined through the dimensionless numbers and are needed in the B.C. definitions. Depending of the problem at hand, a different amount of them will be necessary, but the Reynolds number and Mach number will always have to be defined. In order to have a proper understanding about the non-dimensionalization procedure used in the code, see (@ref mmNonDim).

Special attention when using the **detailed chemistry** equation set has to be taken. Since this equation system has not been non-dimensionalized the infinity values are specified directly and not through the dimensionless numbers.

<table>
<tr><th>Property                         </th><th>Type  </th> <th>Description</th></tr>
<tr><td>Re          </td><td>Float   </td> <td>Reynolds number</td></tr>
<tr><td>Ma          </td><td>Float   </td> <td>Mach number</td></tr>
<tr><td> angle          </td><td> Array of Floats   </td> <td> Angle of the incoming infinity velocity. For example a value [90, 0] would equal to an infinity velocity defined in **Y+** direction. </td></tr>
</table>

## Initial Conditions

Initial conditions describe the starting values of all your variables at the first timestep of your simulation. The closer they are to reality, the faster the convergence and stability of your simulation will be. The initial conditions are implemented in the code, and tend to describe the distribution for a given specific problem. It is assigned by a numerical value to the property `initialCondition` in fvcartesiansolverxd.cpp.  

If the case to simulate has not been performed to date, it may be necessary to specify a new initial condition in the code.

## Numerical settings

#### Surface-value reconstruction

The surface value reconstruction is given by the property `surfaceValueReconstruction`. As the name implies, this property controls what type of algorithm is used to reconstruct the primitive variables at the surfaces from the values at the cell centroids.

Depending on your problem at hand, for example for supersonic simulations where shock-capturing is required, it may be important to switch from the default `"HOCD"`. In this situation, a limited approach may be beneficial to avoid spurious oscillations across the shock.

<table>
<tr><th>Property                         </th><th>Type  </th> <th>Description</th></tr>
<tr><td> surfaceValueReconstruction         </td><td> String </td> <td> Selects the reconstruction method. \n \n
      Possible values are:
      <ul>
      <li> HOCD: unlimited high-order surface data reconstruction. </li>
      <li> HOCD_LIMITED: limited high-order surface data reconstruction. </li>
      <li> LOCD: unlimited low-order (pressure)/high-order surface data reconstruction. </li>
      </ul> </td></tr>
</table>

#### RK-Method
A five-step Runge-Kutta Scheme of second-order accuracy is commonly applied to the temporal integration of the conservative variables. The number of coefficients of the RK method have to be equal to the number of steps specified. An example of the properties is given in the following, where the chosen RK order is set to two, a total of five iteration steps is conducted and the coefficients are chosen in such a way to maximize numerical stability.

````
noRKSteps = 5
rkalpha-step = [0.25, 0.16666666666, 0.375, 0.5, 1.0]
rungeKuttaOrder = 2

````

<table>
<tr><th>Property                         </th><th>Type  </th> <th>Description</th></tr>
<tr><td>noRKSteps          </td><td>Int   </td> <td>Number of steps in the Runge-Kutta time-stepping method.</td></tr>
<tr><td>rkalpha-step          </td><td>   </td> <td> Coefficients of the Runge-Kutta time-stepping method. Possible values are:
    <ul>
    <li>Floating point numbers (as many as Runge-Kutta steps).</li>
    </ul></td></tr>
<tr><td>rungeKuttaOrder          </td><td> Int   </td> <td> Defines the runge kutta method (order). Possible values are:
    <ul>
    <li>2 - second order</li>
    <li>3 - third order</li>
    </ul> </td></tr>
</table>

#### Time step size and CFL number
The chosen time step has a great influence on the numerical stability of the employed methods. It is chosen by assigning a non-negative, non-zero value to the CFL-number `cfl` (Courant–Friedrichs–Lewy number). This number governs the necessary condition for convergence when solving partial differential equations with explicit time integration schemes. Usual values in the theoretical sense range from 0 to 1 for a stable computation, however depending on the method used to compute the time step from the CFL number bigger numbers than 1 can also produce stable solutions. It is encouraged to try different values and find the highest value that still produces stable solutions, since this will speed up the convergence of the simulation while not negatively affecting the accuracy of the simulation if performed properly.

The exact method used to compute the time step from the given CFL number can be defined by the property `timeStepMethod` by assigning an integer value to it. For exact details about how each of the methods works, refer to the code found in fvcartesiansolverxd.cpp. 

<table>
<tr><th>Property                         </th><th>Type  </th> <th>Description</th></tr>
<tr><th>cfl                         </th><th>Float  </th> <th> CFL number</th></tr>
<tr><td> timeStepComputationInterval         </td><td> Int </td> <td> Specifies on which interval the time-step will be recomputed.
    <ul>
    <li> -1 -> default (when restart): never - read from restart file, requires m_restart == true. </li>
    <li> 0  -> default (when no restart): once at the beginning </li>
    <li> 1  -> every time-step </li>
    <li> N  -> every N time-steps </li> </td></tr>
<tr><td> timeStepMethod                        </td><td> Int  </td> <td> Time step computation method. 
    <ul>
    <li> 1 -> computes the time step based on the CFL condition of the Euler equations </li>
    \f$ \Delta t = \min \left[ C * \dfrac{\Delta x}{|u_i| + a} \right], \quad i = 0..d \f$
    <li> 6 -> time step method based on the cell length of the GCells. For level set only. </li>
    <li> 1000 -> time step method for enthalpy solver</li>
    <li> 17511 -> use flame based time steps but write the time step of case 1 in the solution files</li>
    <li> 100 -> Same as 1, but uses the infinity values and fulfills additionally the viscous cfl number </li>
    </ul>
    </td></tr>
</table>

#### Numerical diffusion
Numerical diffusion can be fine tuned through the property `globalUpwindCoefficient`, which sets the upwind coefficient for the pressure splitting 

This property is normally used with values close to zero (around 0.01), meaning vanishing numerical diffusion and could lead to non-physical, high-frequent pressure oscillations for certain simulations. A case by case study may be necessary.

#### Viscous flux computation
The numerical method for the computation of the viscous flux can impact numerical stability of the problem. By default, the five point discretization scheme is used. 

<table>
<tr><th>Property                         </th><th>Type  </th> <th>Description</th></tr>
<tr><td> viscousFluxScheme       </td><td> String </td> <td> Scheme for the calculation of the viscous flux\n
    Possible values are:
    <ul>
      <li>THREE_POINT </li>
      <li>FIVE_POINT </li>
      <li>FIVE_POINT_STABILIZED </li>
    </ul> </td></tr>
</table>

## Memory settings

These settings control the memory allocation of the solver for each computational rank. These values can heavily impact code performance if not chosen appropiately. Trial and error may be required for fine-tuning memory allocation for a given case. A detailed explanation on these settings can be found in (@ugGeneralApplication).

## Boundary condition

Boundary conditions govern the behaviour of the solution in the boundary regions. A distinction is made in the solver between:  

* **regular boundary conditions** : are applied directly at the boundary and computed with the help of external ghost cells.  

* **cut-off boundary conditions** : are applied to surfaces that have been previously "cut-off" the grid. 

Common boundary conditions used in CFD, such as inlet, outlet, wall and symetry are present in the code. More exotic B.C., such as the Navier-Stokes Characteristic Boundary Condition (NSCBC) can also be used. This type of B.C. allows outgoing pressure waves to pass through the boundary without reflecting them, and are useful for acoustic applications or pressure-sensitive simulations.  

A list of the available boundary conditions is generated here:
@subpage fvBcOverview.

The regular boundary conditions are specified in the geometry file by assigning a numerical value (corresponding to the dessired B.C.) to the desired surface. Cut-off B.C. are specified in the properties file (after cut-off treatment) by the following properties:

#### Cut-Off treatment  

<table>
<tr><th>Property                         </th><th>Type  </th> <th>Description</th></tr>
<tr><td>cutOff          </td><td>Boolean   </td> <td>Flag to switch the cut-off treatment</td></tr>
<tr><td>cutOffMethod          </td><td>String   </td> <td>Method to define the cut-off. Accepted values are: 
    * *B*: Box (keeps everything inside the box)
    * *iB*: inverse Box (keeps everything ouside the box)
    * *P*: Plane 
    * *C*: Cylinder</td></tr>
<tr><td>cutOffCoordinates          </td><td>Array of Floats   </td> <td>Coordinates to define the cut-off positions. Definition depends on the chosen method.</td></tr>
</table>

#### Cut-Off B.C.  

<table>
<tr><th>Property                         </th><th>Type  </th> <th>Description</th></tr>
<tr><td>createBoundaryAtCutoff          </td><td>Boolean   </td> <td>Switch to turn the creation of cut-off B.C. on and off</td></tr>
<tr><td>cutOffDirections          </td><td>Array of Ints   </td> <td> Defines the cartesian direction of the cut-off surfaces away from the domain. This means in which direction the cut-off can be found. Following convention is held:
    * **0: X-**
    * **1: X+**
    * **2: Y-**
    * **3: Y+**
    * **4: Z-**
    * **5: Z+**</td></tr>
<tr><td>cutOffBndryIds          </td><td>Array of Ints   </td> <td>The cut-off B.C. associated with the directions given in `cutOffDirections`.</td></tr>
</table>

## Sponge layer
A sponge layer is an artificial modification of the RHS of the discretized conservation equation, such that the difference between the current value of a variable and a given value is dampened. There exists the possibility to define sponge layers around existing, already defined B.C. or to define them on the general boundaries of the domain. Depending on the used sponge layer type, one or more of the variables are dampened. 

For setting the boundaries on the general boundaries of the domain, a `spongeLayerThickness` > 0 has to be defined. The `spongeLayerThickness` is given as a factor of \f$L_{ref}\f$ and controls the thickness of the sponge layer in which the sponge layer forcing is applied. The `spongeLayerType` 
gives the type of sponge chosen for the forcing. A comprehensive list of sponge layer types can be found @subpage fvSponge "here". 

The property `sigmaSponge` controls the amplitude of the forcing term. Bigger numbers should lead a stronger dampening of the difference. The `spongeFactor` is an array and consists of a total of 2 * nDim floating point numbers. It controls on which domain boundaries a sponge layer is active. Each entry corresponds to the respecive direction (0: -x, 1: +x, 2: -y, ...) and controls the specific sponge layer thickness on this domain boundary. If a factor is zero, no sponge layer is generated in this direction.

```  
spongeLayerThickness = 1.0  
spongeLayerType = 91900  
sigmaSponge = 15.0  
spongeFactor = [0.5, 0.5, 2.0, 1.0, 0.5, 0.5] 
```

## RANS
A RANS simulation is started by setting  

`fvSystemEquations = FV_SYSEQN_RANS`

RANS-specific settings:

<table>
<tr><th>Property                         </th><th>Type  </th> <th>Description</th></tr>
<tr><td>ransMethod          </td><td> String   </td> <td> Chooses the RANS equations to be solved. </td></tr>
<tr><td>fullRANS          </td><td> Boolean   </td> <td> Triggers a full RANS simulation.</td></tr>
<tr><td>noRansEquations          </td><td>Integer   </td> <td> </td></tr>
</table>

A mathematical description of the turbulence modeling is given [here](@ref TurbulenceModel.md). Currently, the Spalart-Allmaras turbulence model is implemented. To use it set ```ransMethod = "RANS_SA"``` and ```noRansEquations = 1```.

## Adaptive mesh refinement
If you want to learn how to enable AMR for you simulation, read [here](@ref ugAMR). The list of all available sensors for the FV solver can be found @subpage sensorsFV "here".
