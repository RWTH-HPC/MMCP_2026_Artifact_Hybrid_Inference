# General application control # {#ugGeneralApplication}

## Available solvers
@maia is a numerical simulation framework comprising the solvers listed below. More information on the individual solvers can be found on the respective page.  
1. @ref ugFV
2. @ref ugFVMB
3. @ref ugLB
4. @ref ugDG
5. @ref ugLPT
6. @ref ugLS
7. @ref ugACA

The solvers can be arbitrarily coupled by a coupling class, in which the transfer of information between the solver interfaces is defined (cf. @ref nmCoupling).
The solvers run on subsets of a shared unstructured Cartesian mesh, which is generated in the @ref ugMeshGeneration. This ensures an efficient coupling of the solvers in multisolver simulations. For single solver @ref ugFV simulations also an option to use a @ref ugSTRCD Grid is available.

## Global properties
Some properties are solver independent and are required for each simulation. They describe, for example, for how many time steps a simulation should run, which solvers are used and in which order they are executed. The execution order of the solvers is handled in @ref dvRecipe.
Below an example property file is shown for a coupled @ref ugLS and @ref ugFVMB simulation. A list of all solver independet properties is given in the table on the bottom of the page.

```
### These are the global properties
nDim = 3
gridGenerator = false
flowSolver = true
timeSteps = 100
restartFile = false
restartInterval = 10
restartInterval = 50

outputDir = 'out/'
solutionOutput = './out/'
gridInputFileName = 'grid.Netcdf'
geometryInputFileName = "geometry.toml"

scratchSize = 30

maxNoCells = 10000
minLevel = 5
maxUniformRefinementLevel = 5
maxRfnmntLevel = 6
maxBoundaryRfnLvl = 6

### The properties below specify which solvers and couplers are used
multiSolverGrid = True
noSolvers = 2
solverType_0 = "MAIA_LEVELSET_SOLVER"
solverType_1 = "MAIA_FV_MB"
noCouplers = 1
couplerType_0 = "COUPLER_LS_FV_MB"
solversToCouple_0 = [0,1]
executionRecipe = "RECIPE_INTRASTEP"
recipeMaxNoSteps = 2
solverOrder_0 = [1,0]
solverOrder_1 = [0,1]
couplerOrder_0 = [1,1]
adaptationOrder = [0,1]
```

## Mesh adaptation
The Cartesian mesh can be adapted according to the solver solutions during or before the solver run. Mesh adaptation means that the mesh resolution can be adjusted during the solver run, which can be useful when, for example, a distinct flow feature, such as a traveling shock wave, is supposed to be sufficiently resolved. More information on the mesh adaptation can be found here @ref ugAMR.

## Dynamic load balancing
Due to the mesh adaptation, the computational workload on the individual processors changes during the simulation. To avoid large load imbalances a @ref ugDLB is required.

## List of properties
<table>
<tr><th>Property              </th><th>Type  </th> <th>Description                                 </th></tr>
<tr><td>nDim                  </td><td>Int   </td> <td>Number of space dimensions. </td></tr>
<tr><td>multiSolverGrid       </td><td>Bool  </td> <td>Set to true if more than one solver is used. </td></tr>
<tr><td>gridGenerator         </td><td>Bool  </td> <td>If true the @ref ugMeshGeneration is executed. Musst be unequal to flowSolver. </td></tr>
<tr><td>flowSolver            </td><td>Bool  </td> <td>If true a solver run is started. Musst be unequal to gridGenerator. </td></tr>
<tr><td>postProcessing        </td><td>Bool  </td> <td>Activates the post-processing. </td></tr>
<tr><td>scratchSize           </td><td>Int   </td> <td>During a simulation the memory is dynamically allocated and cleard if required by the solution algorithms. Scratch size defines the buffer size of memory reserved for this dynamic memory allocation. </td></tr>
<tr><td>timeSteps             </td><td>Int   </td> <td>The number of time steps that will be executed. </td></tr>
<tr><td>testcaseDir           </td><td>String</td> <td>Path to simulation directory (Best leaf at "./"). </td></tr>
<tr><td>restartDir            </td><td>String</td> <td>Directory from which restart files are read in case of a simulation restart. Use relative path from testcaseDir. If empty restartDir is set to outputDir. </td></tr>
<tr><td>inputDir              </td><td>String</td> <td>Input directory for ACA solver and location of geometry file (geometryInputFileName). Use relative path from testcaseDir. </td></tr>
<tr><td>outputDir             </td><td>String</td> <td>Directory to which restart files are written. Use relative pwth from testcaseDir. </td></tr>
<tr><td>solutionOutput        </td><td>String</td> <td>Directory to which solver specific solution files are written. </td></tr>
<tr><td>gridInputFileName     </td><td>String</td> <td>Name of grid file to be used for the simulation. </td></tr>
<tr><td>geometryInputFileName </td><td>String</td> <td>Name of geometry file to be used for the simulation. </td></tr>
<tr><td>restartFile           </td><td>Bool  </td> <td>Specifies whether a simulations is newly initialized or restarted from specific restart file. </td></tr>
<tr><td>restartTimeStep       </td><td>Int   </td> <td>Time step from which simulation is restarted (requires restartFile=true). </td></tr>
<tr><td>restartInterval       </td><td>Int   </td> <td>Interval in which restart files are written. </td></tr>
<tr><td>solutionInterval       </td><td>Int   </td> <td>Interval in which solver solution files are written. </td></tr>
<tr><td>useNonSpecificRestartFile</td><td>Bool</td><td>If true, previous restart files will be overwritten and current time step is read from restart file. Property 'restartTimeStep' will be ignored. </td></tr>
<tr><td>noSolvers             </td><td>Int   </td> <td>The number of solvers used in the simulation. </td></tr>
<tr><td>solvertype            </td><td>String</td> <td>The type of the solver to be created. The solvers grid type is automatically determined. </td></tr>
<tr><td>noHaloLayers          </td><td>Int   </td> <td>The number of halo layers for parallel simulations or periodic domains. The number of halo layers between solvers can differ. </td></tr>
<tr><td>noCouplers            </td><td>Int   </td> <td>The number of couplers used for the coupling of the solver. </td></tr>
<tr><td>couplerType_          </td><td>String</td> <td>The type of coupler that is used. The coupler type depends on the solvers that are coupled. </td></tr>
<tr><td>executionRecipe       </td><td>String</td> <td>Sets the type of the execution @ref dvRecipe. </td></tr>
<tr><td>recipeMaxNoSteps      </tb><td>Int   </td> <td>This property specifies the MAIAExecutionRecipe::m_maxNoSteps of sub steps performed in the execution @ref dvRecipe.</td></tr>
<tr><td>solverOrder_*         </td><td>Int   </td> <td>This property if a solver is active during a specific sub step of the execution @ref dvRecipe. solverOrder needs to be set for each solver individually, where * is replaced by the solver id. This property is given as an array, e.g., for MAIAExecutionRecipe::m_maxNoSteps = 2: [1,0] if solver is active in first step and inactive in second step. </td></tr>
<tr><td>couplerType_*         </td><td>String</td> <td>Sets the type of the couplers. The couplerType depends on the solvers used. This property is set for each coupler indivudally by replacing * with the coupler id. </td></tr>
<tr><td>couplerOrder_*        </td><td>Int   </td> <td>This property if a solver is active during a specific sub step of the execution @ref dvRecipe. solverOrder needs to be set for each solver individually, where * is replaced by the solver number. This property is given as an array, e.g., for MAIAExecutionRecipe::m_maxNoSteps = 2: [1,0] if solver is active in first step and inactive in second step. </td></tr>
<tr><td>solversToCouple_*     </td><td>Int   </td> <td>Detremines which solver are coupled by the specific coupler. Replace * with the respective coupler id. Each coupler couples to solvers. E.g., write [1,2] to couple solvers 1 and 2. </td></tr>
<tr><td>maxNoCells            </td><td>Int   </td> <td>For each grid cell handled by a single processor, memory allocated. This property limits the maximum number of cells per processor and thereby influences the memory foot print. </td></tr>
<tr><td>minLevel              </td><td>Int   </td> <td>The grid level on which the partitioning of the grid is performed if no partition level shift exists. The minLevel is set during the @ref ugMeshGeneration. Maybe you can find some information on partition level shifts in section @ref ugHPC or @ref nmParallelization. </td></tr>
<tr><td>maxUniformRefinementLevel</td><td>Int</td> <td>The maximum grid level to which each minLevel cell is refined. See also @ref ugMeshGeneration. </td></tr>
<tr><td>maxBoundaryRfnLvl     </td><td>Int   </td> <td>The maximum grid level to which cells at the speciefied domain boundaries are refined. </td></tr>
<tr><td>maxRfnmntLvl          </td><td>Int   </td> <td>The maximum grid level allowed. </td></tr>
<tr><td>adaptation            </td><td>Bool  </td> <td>If true, the @ref ugAMR  of the Cartesian mesh is allowed. </td></tr>
<tr><td>initialAdaptation     </td><td>Bool  </td> <td>If a new simulation is initialized (restartFile=false), a mesh adaptation based on the solver intital conditions can be forced. </td></tr>
<tr><td>adaptationOrder       </td><td>Int   </td> <td>Specifies in which sub steps of the execution @ref dvRecipe a mesh adaptation is allowed. Specify as array, e.g., for MAIAExecutionRecipe::m_maxNoSteps = 2: [1,0] if adapation is only allowed after the first sub step. </td></tr>
<tr><td>balance               </td><td>Bool  </td> <td>If true, the @ref ugDLB is activated. </td></tr>
<tr><td>dualTimeStepping      </td><td>Bool  </td> <td>The #ref ugFV solver featrues a dual time stepping method. (Untested, use at own risk.) </td></tr>
<tr><td>splitMpiComm          </td><td>Bool  </td> <td>See @ref ugFV and @ref ugDG for more information. If unsure set splitMpiComm=false. </td></tr>
<tr><td>noPostProcessing      </td><td>Int   </td> <td>The number of post-processing solvers. For more information on post-processing see @ref ugPostprocessing. </td></tr>
</table>

