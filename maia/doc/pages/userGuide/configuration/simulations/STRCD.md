# Structured Finite Volume (FVSTRCTRD) # {#ugSTRCD}

[TOC]

## Overview

m-AIA is divided into two distinct parts, one uses a
hierarchical Cartesian cut-cell approach with automatic grid
generation techniques.  The other part is the structured block
environment for solving the Navier-Stokes equation (NSE) on
curvilinear boundary fitted grid. The @ref nmStructuredGrid needs to
be provided beforehand as the solver is not able to create them. This
page gives introduction on how to configure your properties files for the
structured finite volume (FVSTRCTRD) solver environment.

At present, the FVSTRCTRD solver can be used for both RANS and LES, also
in coupled simulation using RANS for some blocks and LES for
others. Simulations can be conducted in 2D and 3D for static and
moving meshes.

## General settings 

Some common properties, described in @ref ugGeneralApplication, are also necessary for the FVSTRUCTURED solver, including:

```
### These are the global properties
nDim = 3
flowSolver = true
timeSteps = 100
restartFile = false

outputDir = 'out/'
solutionOutput = './out/'
gridInputFileName = 'grid.hdf5'

scratchSize = 30
maxNoCells = 10000

### The properties below specify which solvers and couplers are used
multiSolver = 1
noSolvers = 1
executionRecipe = "RECIPE_BASE"
solvertype.default = "MAIA_UNIFIED"
solvertype.0 = "MAIA_STRUCTURED"

### The following are not required for FVSTRUCTURED but the general framework still reads them
gridGenerator = false
maxRfnmntLvl = 8
```

After this, you can choose LES or RANS

### DNS/LES

DNS/Implicit LES is the default case, you do not need to set any properties to turn it on.

To deactivate the viscous flux computation you can use:
~~~
govEqs = "EULER"
~~~

### RANS

To activate RANS, it is necessary to make the following settings

~~~
fullRANS = true
~~~

Choose the desired RANS method:
~~~
ransMethod = "RANS_SA_DV"
~~~

## Initial condition
The initial condition is specified via a number

~~~
initialCondition = 1233
~~~

where the specific IC is found in the following table.

<table>
<tr> <th> initialCondition   </th> <th> Description </th>   <th> Applicable for RANS </th></tr>  
<tr> <td> 0                    </td> <td> Parallel inflow field. </td>   <td> Yes  </td> </tr>
<tr> <td> 2                   </td> <td> Shear flow       </td>   <td> No              </td> </tr>
<tr> <td> 11                   </td> <td> Point source in the middle      </td>   <td> No              </td> </tr>
<tr> <td> 43/333               </td> <td> Parallel inflow field with pressure peak in the middle of the domain. </td>   <td> No              </td> </tr>
<tr> <td> 314                  </td> <td> Stagnating flow field \f$ u = v = w = 0 \f$, \f$\rho = \rho_\infty\f$, \f$p = p_\infty\f$        </td>   <td> No             </td> </tr>
<tr> <td> 315                  </td> <td> Poiseuille flow        </td>   <td> No              </td> </tr>
<tr> <td> 101                  </td> <td> Taylor-Green Vortex   \f$\rho = \rho_\infty\f$       </td>   <td> No              </td> </tr>
<tr> <td> 1234                 </td> <td> Laminar channel flow          </td>   <td> No              </td> </tr>
<tr> <td> 1233                 </td> <td> Turbulent channel flow with perturbations          </td>   <td> No             </td> </tr>
<tr> <td> 1236                 </td> <td> Pipe flow with perturbations         </td>   <td> No              </td> </tr>
<tr> <td> 79092                </td> <td> Approximate mean turbulent boundary layer profile         </td>   <td> Yes              </td> </tr>
<tr> <td> 79091                </td> <td> Turbulent plate with log-law          </td>   <td> Yes              </td> </tr>
<tr> <td> 111/112/113/44          </td> <td> Three-dimensional non-reflective boundary conditions  </td>   <td> No              </td> </tr>
<tr> <td> 777                </td> <td> Falkner Skan Cooke initial condition (incompressible)         </td>   <td> No              </td> </tr>
<tr> <td> 999                </td> <td> Blasius laminar boundary layer       </td>   <td> No              </td> </tr>
</table>

## Boundary condition
At the boundaries of the domain suitable boundary conditions need to be applied. These are specified in the grid file by assigning numbers to each window for each block:
<table>
<tr><th> BC </th> <th>Description</th></tr>

<tr><td> 1000 </td> <td> **No-slip wall**
* \f$ u_{wall} = v_{wall} = w_{wall} = 0\f$
* \f$\left.\frac{\partial \rho}{\partial n}\right|_{wall} = 0 \f$
* \f$\left.\frac{\partial p}{\partial n}\right|_{wall} = 0 \f$

When moving grid is switched on this will become 1004.
</td> </tr>

<tr><td> 1003 </td> <td> **Isothermal no-slip wall**
* \f$ u_{wall} = v_{wall} = w_{wall} = 0\f$
* \f$ T_{wall}  = \alpha  T_\infty \f$
* \f$\left.\frac{\partial p}{\partial n}\right|_{wall} = 0 \f$

When moving grid is switched on this will become 1006. For this BC, **isothermalWallTemperature** in the properties files sets the value of \f$ \alpha \f$ which controls the wall temperature as a ratio of the freestream temperature.
</td> </tr>

<tr><td> 1004 </td> <td> **Moving no-slip wall**
* \f$ u_{wall} = u_{mesh}(t)\f$
* \f$ v_{wall}  = v_{mesh}(t)\f$
* \f$ w_{wall} = w_{mesh}(t)\f$
* \f$ \left.\frac{\partial \rho}{\partial n}\right|_{wall} = 0 \f$
* \f$ \left.\frac{\partial p}{\partial n}\right|_{wall} = 0 \f$.

That is, the velocity of the boundary grid points, according to the prescribed boundary motion, are applied to the fluid. </td> </tr>

<tr><td> 1006 </td> <td> **Moving isothermal no-slip wall**
* \f$ u_{wall} = u_{mesh}(t)\f$
* \f$ v_{wall}  = v_{mesh}(t)\f$
* \f$ w_{wall} = w_{mesh}(t)\f$
* \f$ T_{wall}  = \alpha * T_\infty \f$
* \f$ \left.\frac{\partial p}{\partial n}\right|_{wall} = 0 \f$.

That is, the velocity of the boundary grid points, according to the prescribed boundary motion, are applied to the fluid. For this BC, **isothermalWallTemperature** in the properties files sets the value of \f$ \alpha \f$ which controls the wall temperature as a ratio of the freestream temperature. </td> </tr>

<tr><td> 2002 </td> <td> **Supersonic inflow**
* \f$ u = u_{\infty}, v = v_{\infty}, w = w_{\infty}\f$
* \f$ \rho = \rho_{\infty}, p = p_{\infty}.\f$
</td> </tr>

<tr><td> 2003 </td> <td> **Simples subsonic in/outflow**

Simple extrapolation and prescription

For inflows, i.e., \f$M < 0\f$, the following conditions are used

* \f$ u = u_\infty\f$
* \f$ v = v_\infty\f$
* \f$ w = w_\infty\f$
* \f$\frac{\partial \rho}{\partial n} = 0 \f$
* \f$\frac{\partial p}{\partial n} = 0 \f$

On outflow boundaries, i.e., where \f$M > 0\f$, the conditions


* \f$\frac{\partial u}{\partial n}  = 0 \f$
* \f$\frac{\partial v}{\partial n}  = 0 \f$
* \f$\frac{\partial w}{\partial n}  = 0 \f$
* \f$\frac{\partial \rho}{\partial n}  = 0 \f$
* \f$ p = p_\infty \f$

</td> </tr>

<tr><td> 2004 </td> <td> **Subsonic in/outflow**

Depending on the flow direction this acts as an in- or outflow, simplified characteristics approach, Whitfield 1984

For inflows, i.e., \f$M < 0\f$, the following conditions are used
\f{align}{
  p|_{BC} &= \frac{1}{2} \left[p|_i + p_\infty + \frac{\rho_{old} a_{old}}{\left| \nabla \vec{\xi}\right|}\left(\vec{\xi}\cdot\left(\vec{u}|_i - \vec{u}_\infty\right)\right)\right]\\
  \rho|_{BC} &= \rho_\infty + \frac{1}{a_{old}^2}\left(p|_{BC} - p_\infty\right)\\
  \vec{u}|_{BC} &= \vec{u}_\infty + \frac{\vec{\xi}}{\left|\nabla \vec{\xi}\right|\rho_{old}a_{old}}\left(p|_{BC} - p_\infty\right)\,,
\f}

On outflow boundaries, i.e., where \f$M > 0\f$, the conditions
\f{align}{
  p|_{BC} &= p_\infty\\
  \rho|_{BC} &= \rho|_i + \frac{1}{a_{old}^2}\left(p|_{BC} - p|_i \right)\\
  \vec{u}|_{BC} &= \vec{u}|_i - \frac{\vec{\xi}}{\left|\nabla \vec{\xi}\right|\rho_{old}a_{old}}\left(p|_{BC} - p|_{i}\right)
\f}

</td> </tr>
<tr><td> 2005 </td> <td> **Supersonic outflow**

* \f$\frac{\partial u}{\partial n} = 0 \f$
* \f$\frac{\partial v}{\partial n} = 0 \f$
* \f$\frac{\partial w}{\partial n} = 0 \f$
* \f$\frac{\partial p}{\partial n} = 0 \f$
* \f$\frac{\partial \rho}{\partial n} = 0 \f$

</td> </tr>
<tr><td> 2500,2501 </td> <td> **Recycling-rescaling**


A combined rescycling-rescaling boundary condition for compressible
wall-bounded flows by [El-Askary et al.][ElAskary2003]. Takes a slice
from the window connected to BC 2501 (needs to be specified in the
grid file, typically at some distance downstream of the inflow
window), rescales it to the desired boundary layer thickness, and
prescribes it to the inflow window with BC 2500 to achieve a turbulent
boundary layer flow.  The property **rescalingBLT** in the properties
file is used to rescale the boundary layer at the inflow which means
the smaller this value is the thinner the boundary layer is. The heart of this method is a
means of estimating the velocity at the inlet plane based on the
solution downstream. The velocity field is extracted from a plane near
the domain exit, which is uesed to rescale, and then it is
reintroduced as a boundary condition at the inlet.

</td> </tr>

<tr><td> 2600,2601 </td> <td> **External profile**

Prescribe a given profile to the window from an external HDF5 file.
</td></tr>

<tr><td> 4401,4403,4405 </td> <td> **Periodic boundary**

Periodic boundary condition for two opposing windows. Use the same odd number for two opposing windows, for multiple periodic directions use different numbers, e.g. 4401 for the first two windows (first periodic direction), 4403 for other two opposing windows (second periodic direction), etc.
</td> </tr>

<tr><td> 6000 </td> <td> **Multiblock connection**

Enables Multiblock communication between user defined blocks in the grid file
</td> </tr>

<tr><td> 7909 </td> <td> **Synthetic turbulence generation (STG)**

A modified version of the reformulated synthetic turbulence generation(RSTG) method by [Roidl et al.][Roidl2013] is used.
Coherent structures are generated by superimposing the velocity fluctuations of virtual eddy cores, which populate  a virtual volume, onto the inflow plane.

The synthetic turbulence generation inflow boundary condition is used by setting the inflow window boundary condition number to **7909** and setting the property

~~~
useSTG = true
~~~

At the initial startup, which is indicated by setting
~~~
stgInitialStartup = true
~~~

The averaged primitive variables \f$ (\overline{u}, \overline{v}, \overline{w}, \overline{p}, \overline{\rho})\f$ and the turbulent viscosity \f$ \overline{\nu_t}\f$ need to be in the first layers of cells near the inflow. This can be achieved, e.g., by interpolating the result from a precursor RANS simulation to the flow field of the current LES case. The boundary layer thickness \f$ \delta_{99} \f$ of the average velocity profile at the inflow needs to be set via

~~~
deltaIn = 7.08
~~~

and the number of eddies is set by


~~~
stgMaxNoEddies = 100
~~~

where the number of eddies should be \f$ n_{eddies} \approx L_z/\delta_{99}*100 \f$. The property

~~~
stgCreateNewEddies = false
~~~

can normally be left to false, only when the generation of a completely new distribution is desired this property can be switched on once.

The properties

~~~
BLT1 = 1.0
BLT2 = 1.0
BLT3 = 1.0
~~~

control the size of the virtual box by specifying a factor in each space dimension by which the given \f$\delta_{99}\f$ is multiplied to compute the size in the \f$x\f$- and \f$y\f$-direction, in the spanwise direction \f$z\f$ the factor is multiplied with the spanwise size of the domain. For instance, using the factors above the virtual box extends one \f$\delta_{99}\f$-thickness in the \f$x\f$-drection, centered around the inflow plane, one \f$\delta_{99}\f$-thickness in the \f$y\f$-direction, starting at the bottom wall, and covers exactly the spanwise domain extent.

The property

~~~
stgExple = 0.5
~~~

controls the exponent in the length scale computation, typically a value of 0.5 delivers reasonable results.

 </td> </tr>
</table>


## Output

In MAIA there are different options to write out your solutions: full solution output, box output, interpolated output, and force field output.

### Solution/Restart output

The normal solution output saves the solutions of the whole flow field, these files can also be used as restart files. That is, in the FVSTRCTRD solver, any solution file can be used as a restart file.  
|property | Definition|
| ------ | ------ |
|outputDir | Directory to which the files will be written |
|solutionInterval | Interval of output in time steps |

To restart the computation use the following properties:

| property | Definition|
| ------ | ------ |
| restartFile | Set true for restart  |
| restartVariablesFileName| Name of the restart file|


### Boxes output
This output can be used to write out 3D subvolumes of the domain. The offsets and sizes of the subvolume are given in computational indices \f$(i,j,k)\f$. An example is found in the *3D_prescribingBC_LES* testcase. The properties are:

|property | Definition|
| ------ | ------ |
|boxOutputInterval | Number of timesteps between writing out the files |
|boxBlock | Number of boxes to be written out |
|boxOffsetK | cell-offset in k-direction |
|boxOffsetJ | cell-offset in j-direction |
|boxOffsetI | cell-offset in i-direction |
|boxSizeK | size in cells in k-direction |
|boxSizeJ | size in cells in j-direction |
|boxSizeI | size in cells in i-direction |
|boxWriteCoordinates | set true to write out cell-center coordinates |

@note
If there is no box output wanted, set the **boxOutputInterval** to 0.

### Forces output

The forces on all surfaces with wall boundary conditions, i.e., any boundary condition \f$ 1000 \geq 1000 < 2000\f$, can be written instantaneously to HDF5 files (full 2D resolution) with name **auxData????.hdf**, where ???? is the time step, or integrated to scalar coefficients to a text file with name **forces.?.dat**, where ? is the number of wall for which the coefficients were computed. That is, if you have three surfaces with a wall boundary condition three text files __forces.0.dat__, __forces.1.dat__, and __forces.2.dat__ are written.

To trigger the 2D output to HDF5 files set the following property:
~~~
computeForce              = true                  
~~~

This will trigger the force calculation, inclusive skin friction and pressure forces. For moving walls the power consumption can also be computed by setting
~~~
computePower              = true
~~~

The force coefficients, which are computed by integrated over the forces over the wall surface, are written as attributes to the HDF5 files. Averaging the coeffient along \f$i\f$, \f$j\f$, or \f$k\f$ can be triggered by
~~~
computeCpLineAverage      = false
cpAveragingDir            = 0
~~~

where __cpAveragingDir__ is set to \f$0\rightarrow i\f$, \f$1\rightarrow j\f$, or \f$2\rightarrow k\f$.

The parameter below will define the output interval of the auxdata HDF5 files.
~~~
forceOutputInterval       = 10000
~~~

To trigger the writeout of the integrated force coefficients to a text file you can set

~~~
forceAsciiOutputInterval  = 50
forceAsciiComputeInterval = 1
~~~

where **forceAsciiComputeInterval** controls how often the coefficients are computed, e.g., 1 for every time step, and **forceAsciiOutputInterval** controls how often the data is written out to the text files. For the values above the forces are computed in every time step, but saved internally, and are written out only every 50 time steps to the text files, this procedure will save computational time.

The region where the forces are computed can be limited by

~~~
auxDataCoordinateLimits   = true
auxDataLimits             = [155.0, 205.0, -30.0, 130.0]
~~~

The first two number of **auxDataLimits** control the lower and upper interval boundary in the \f$x\f$-direction and the latter two control the lower and upper interval boundary in the \f$z\f$-direction.

### Interpolated points output

This output interpolates the cell-centered variables to arbitrary points which are located equidistantly on given lines. Currently, one or two directions can be specified, i.e., 1D and 2D distribution of points can be realized. An example of this can be found in the **3D_cube_periodic_LES** testcase. The parameters are:

|property | Definition|
| ------ | ------ |
|intpPointsOutputInterval | Number of timesteps between writing out the files |
|intpPointsStartX | x-Coordinate of the starting point of the line | 
|intpPointsStartY | y-Coordinate of the starting point of the line | 
|intpPointsStartZ | z-Coordinate of the starting point of the line | 
|intpPointsDeltaX | delta in x-Direction between two points on the line |
|intpPointsDeltaY | delta in x-Direction between two points on the line |
|intpPointsDeltaZ  | delta in x-Direction between two points on the line |
|intpPointsDeltaX2D | delta in x-Direction between two points on the line for the second direction|
|intpPointsDeltaY2D | delta in x-Direction between two points on the line for the second direction|
|intpPointsDeltaZ2D  | delta in x-Direction between two points on the line for the second direction|
|intpPointsNoPoints | delta in x-Direction between two points on the line |
|intpPointsNoPoints2D | delta in x-Direction between two points on the line for the second direction |

## Moving grid methods

To activate moving grids set
~~~
movingGrid = true
~~~

The following list of moving grid methods uses hard-coded parameters, which cannot be set in the property file:

<table>
<tr><th>  gridMovingMethod  </th> <th> Description </th> </tr>
<tr><td> 1      </td> <td> Oscillation in \f$i\f$-direction, generating a Stokes layer </td> </tr>
<tr><td> 2      </td> <td> Channel with oscillating indentation      </td> </tr>
<tr><td> 3      </td> <td> Piston moving in \f$i\f$-direction      </td> </tr>
<tr><td> 4      </td> <td> Inner grid movement. The grid points in the center of the domain are actuated in the \f$x\f$- and \f$y\f$-direction with a gradual weakening of the movement towards the boundaries. The mesh at the boundaries is static. </td> </tr>
</table>

In the following table a number of more complex moving grid functions is described which can be controlled via the property file.

### Traveling wave

Spanwise traveling transversal surface wave with parameters defined in inner units (__gridMovingMethod = 9__) or outer units (__gridMovingMethod = 10__)
<table>
<tr><th> properties       </th> <th> Description </th> </tr>
<tr><td> gridMovingMethod                 </td> <td> 9/10 </td> </tr>
<tr><td> waveLengthPlus/waveLength        </td> <td> The wavelength of traveling wave with inner scale/outer scale      </td> </tr>
<tr><td> waveAmplitudePlus/waveAmplitude  </td> <td> The amplitude of traveling wave with inner scale/outer scale      </td> </tr>
<tr><td> waveTimePlus/waveTime      </td> <td> The wavetime of traveling wave with inner scale/outer scale      </td> </tr>
<tr><td> waveBeginTransition     </td> <td> The start of wave spatial transition along i direction, from flat surface to wave surface </td> </tr>
<tr><td> waveEndTransition     </td> <td> The end of wave spatial transition along i direction, from flat surface to wave surface      </td> </tr>
<tr><td> waveOutBeginTransition     </td> <td> The start of wave spatial transition along i direction, from wave surface to flat surface      </td> </tr>
<tr><td> waveOutEndTransition     </td> <td> The end of wave spatial transition along i direction, from wave surface to flat surface      </td> </tr>
<tr><td> waveAngle     </td> <td> The angle of wave     </td> </tr>
<tr><td> waveYBeginTransition     </td> <td> The start of wave spatial transition along j direction, from wave to flat      </td> </tr>
<tr><td> waveYEndTransition     </td> <td> The end of wave spatial transition along j direction, from wave to flat       </td> </tr>
<tr><td> waveTemporalTransition     </td> <td> The time setting for wave temporal trasition  </td> </tr>
</table>


### Oscillating cylinder

This controls the up- and down oscillation of a cylinder, i.e., a circular mesh around a round cylinder is required.

<table>
<tr><th> properties       </th> <th> Description </th> </tr>
<tr><td> gridMovingMethod      </td> <td> 14  </td> </tr>
<tr><td> oscAmplitude      </td> <td> The amplitude of oscillation      </td> </tr>
<tr><td> oscSr      </td> <td> ???      </td> </tr>
<tr><td> oscFreqFactor      </td> <td> The frequency factor used to control the oscillating frequency      </td> </tr>
</table>

## Sponge

To reduce the numerial reflections from boundaries a damping sponge
layer zone can be introduced to absorb outgoing waves.

In the FVSTRCTRD solver the sponge is realized by an additional local
damping term, the strength and type of the damping as well as the
region where the sponge is applied can be controlled via properties.

Activate the sponge by
~~~
useSponge                = true
~~~

The location where sponges are applied is controlled by

~~~
readSpongeFromBC         = false
spongeBndryCndIds        = [7909]
spongeWindowIds          = [1, 2, 4]
~~~

If **readSpongeFromBC** is true the sponge will be added inwards from
all windows which have a BC in the list
**spongeBndryCndsIds**. Otherwise, if **readSpongeFromBC** is false,
the list **spongeWindowIds** will be used to apply the sponge inwards
from these specific window ids (as given in the grid).

The type and strength of the sponge are controlled with the following properties

<table>
<tr><th> properties       </th> <th> Description </th> </tr>
<tr><td> spongeLayerType = 6 </td> <td>

Type of sponge layer:

* 1: sponge to \f$\rho_\infty\f$ and \f$(\rho E)_\infty\f$
* 2: sponge to \f$\overline{\rho}\f$ and \f$\overline{\rho E}\f$ from the fields **spongeRho** and **spongeRhoE** from the restart file
* 3: sponge to \f$\rho_\infty\f$ and \f$p_\infty\f$
* 4: sponge to \f$\overline{\rho}\f$ from the field **spongeRho** from the restart file and \f$ p_\infty\f$
* 5: sponge to \f$p\f$ from the Falkner-Skan Cooke pressure function
* 6: sponge to \f$p_\infty\f$


</td></tr>
<tr><td> spongeLayerThickness = [50.0, 30.0, 20.0] </td> <td> Thickness of the sponge layer normal to the window </td></tr>
<tr><td> sigmaSponge = [1.0, 1.0, 1.0] </td> <td> Sponge strength       </td></tr> 
<tr><td> betaSponge = [2.0, 2.0, 2.0]  </td> <td>  Controls the decay exponent of the sponge strength inwards from the window.  Linear: spongeBeta = 1.  Quadratic: spongeBeta = 2.    </td></tr>
<tr><td> computeSpongeFactor = false  </td> <td> Disables the sponge strength computation. If set to false, the sponge values will be read from the restart file
 which may be much faster due to the slow sponge computation. That is, the sponge needs to be computed once (first solver run, set computeSpongeFactor = true) and is then saved to the restart file. Following simulation runs (set computeSpongeFactor = false) the sponge factor is read from the restart file.</td></tr>
</table>


## Postprocessing

The FVSTRCTRD solver uses a standalone postprocessing part, independent from the other Cartesian-grid solvers.

To activate postprocessing ste
~~~
postprocessing = true
~~~

Averaging of the flow variables can be done while the simulation is running (in solve) without the need to write out snapshots. Otherwise, a set of snapshots written to disk can be read in successively by the solver and averaged. Furthermore, other postprocessing options are available.

The general postprocessing operations which are desired can be controlled by specifiying 
~~~
postprocessingOps = ["PP_AVERAGE_IN", "PP_TAUW_PRE"]
~~~

In the following table are the postprocessing operations which are currently available. Not all operations can be used together with all other operations.

<table>
<tr><th> postprocessingOps </th> <th> Description </th> </tr>
<tr><td> PP_AVERAGE_PRE </td> <td> Averaging of existing solution files before the simulation start         </td></tr>
<tr><td> PP_AVERAGE_POST           </td> <td>     Averaging of existing solution files after the simulation ends     </td></tr>
<tr><td> PP_AVERAGE_IN          </td> <td>  Averaging of the flow variables during the simulation (in solve)       </td></tr>
<tr><td> PP_TAUW_PRE            </td> <td> Compute the average skin friction from the primitive variables </td></tr>
<tr><td> PP_LOAD_AVERAGED_SOLUTION_PRE    </td> <td> Loads the averaged variables from the Mean*.hdf5 file into the solver and into the primitive variables </td></tr>
<tr><td> PP_SUBTRACT_MEAN   </td> <td> Subtracts the mean from the currently loaded solution file, mean file needs to be available  </td></tr>
<tr><td> PP_SUBTRACT_PERIODIC_FLUCTUATIONS  </td> <td> Special option for traveling wave actuation. The periodic fluctuations (from a triple decomposition) are subtracted from the currently loaded solution file </td></tr>
<tr><td> PP_WRITE_GRADIENTS  </td> <td> Computes the spatial gradients of the primitive variables and writes them to an HDF5 file </td></tr>
<tr><td> PP_DECOMPOSE_CF           </td> <td> Perform an FIK decomposition  </td></tr>
<tr><td> PP_COMPUTE_PRODUCTION_TERMS_PRE    </td> <td> Compute production terms  </td></tr>
<tr><td> PP_COMPUTE_DISSIPATION_TERMS_PRE   </td> <td> Compute dissipation terms from a set of solution files and an existing mean file </td></tr>
</table>


General averaging is controlled by the properties

<table>
<tr><th> Property </th> <th> Description </th> </tr>
<tr><td>  pp_averageStartTimestep           </td> <td> Timestep at which averaging is started         </td></tr>
<tr><td>  pp_averageStopTimestep            </td> <td> Timestep at which averaging is ended         </td></tr>
<tr><td>  pp_averageInterval                </td> <td> Interval at which probes of the flow variables are added onto the average   </td></tr>
<tr><td>  pp_fileName           </td> <td> Input file name from which a mean file is read if the postprocessing operation requires it </td></tr>
</table> 

## Interpolation

The flow data for grids which do not match the current grid can be interpolated to the current grid with a built-in method. The interpolation procedure is performed at the initial startup. The interpolation is second-order accurate if a cell center on the current grid is surrounded by 8 cell centers of the donor grid. Then a trilinear interpolation, modified for curvilinear grids, is performed. If no surrounding donor points can be found the fallback procedure resorts to a simple nearest neighbor interpolation. Therefore, for a high quality of the interpolation the goal is to have a donor grid which is at least as big as the current grid, such that all points on the current grid are covered.

It is also possible to interpolate from 2D grids to 3D grids.

Here are properties to be set:
<table>
<tr><th> Property </th> <th> Description </th> </tr>
<tr><td>  restartInterpolation      = true </td> <td>  Turn on this function       </td></tr>
<tr><td>  donorGridName             = "donorGrid.hdf5" </td> <td> The donor grid file         </td></tr>
<tr><td>  donorVars                 = "donorVars.hdf5"</td> <td>  The donor solution file   </td></tr>
<tr><td>  donorScale                = [1.0,1.0,1.0] </td> <td>  Scaling factor of the donor grid  </td></tr>
<tr><td>  donorTranslation          = [0.0,0.0,0.0] </td> <td>  Translation of the donor grid  </td></tr>
</table>

## Zonal computations
High-resolution LES of a turbulent flow at high Reynolds numbers consumes high computational power and time due to the necessary fine grid resolution. 
However, a turbulent flow problem in non-equilibrium state occurs only in a specific area. Applying LES throughout the domain is therefore not an efficient way to solve the turbulent flow.
Therefore, it would be elegant to perform a RANS simulation for noninteresting areas in the flow field, while using the LES simultaneously in the important regions where a higher grid resolution is required.
The FVSTRCTRD solver offers fully coupled zonal RANS-LES simulations using multiple RANS and LES zones. The grid should contain at least two blocks which should overlap by a certain amount. In the property file the type of closure model (RANS or LES) is then specified for each block. The transition from RANS to LES in the near-wall region is typically done by using the STG boundary condition.

Here are relevant properties:
<table>
<tr><th> Property </th> <th> Description </th> </tr>
<tr><td> zonal                      = true  </td> <td>  To turn on zonal computations       </td></tr>
<tr><td> fullRANS                   = false</td> <td>   If you use the zonal approach, fullRANS needs to be false       </td></tr>
<tr><td> noRansZones                = 1  </td> <td> Number of blocks with RANS        </td></tr>
<tr><td> ransZone                   = [0]</td> <td> Block ids (in the grid file) which should be RANS  </td></tr>
<tr><td> ransMethod                 = "RANS_SA_DV"    </td> <td> RANS model to be used  </td></tr>
<tr><td> zonalExchangeInterval      = 50    </td> <td> Timesteps Interval for exchanging data between RANS and LES zones   </td></tr>
<tr><td> zonalAvgFactor             = 128  </td> <td> Averaging factor for the moving average in the LES zones </td></tr>
</table>

## Test cases

For STRCTRDFV solver, there are complete **workshops@ref ...** to start with. Furthermore, some ready-to-use testcases are found in the STRUCTURED_FV subfolder of the testcase repository:
~~~
cd ~/scratch
svn co http://svn/maia/testcases/STRUCTURED_FV/
~~~


## References

* W. El-Aksary, W. Schröder, M. Meinke. LES of compressible wall-bounded flows. AIAA 2003-3554 (2003. [10.2514/6.2003-3554][ElAskary2003]
* B. Roidl, M. Meinke, W. Schröder. A reformulated synthetic turbulence generation method for a zonal RANS-LES method and its application to zero-pressure gradient boundary layers. Int. J. of Heat Fluid Flow, 44 (2013) 28-40. [https://doi.org/10.1016/j.ijheatfluidflow.2013.03.017][Roidl2013]

[ElAskary2003]: http://dx.doi.org/10.2514/6.2003-3554
[Roidl2013]: https://doi.org/10.1016/j.ijheatfluidflow.2013.03.017
