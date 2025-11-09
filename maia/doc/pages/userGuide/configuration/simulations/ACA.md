# Acoustic analogy (ACA) # {#ugACA}

When using the ACA solver the workflow usually consists of three steps
1. [Generating input](@ref ugACA_generatingInput) describing temporal series on
a FWH permeable surface
  1. Generating a surface grid and gathering data during a simulation run such
  as the [FV](@ref ugFV), [LB](@ref ugLB), or [DG](@ref ugDG) solver on this
  surface grid
  2. By describing the source term via an analytical function such as Mono-, Di-,
  or Quadropole by the ACA solver itself
2. [Far-Field prediction](@ref ugACA_farFieldPrediction) for specified observer
locations using the ACA solution methods
3. Optionally [post-process](@ref ugACA_postProcessing) output at observer
location which is dumped in time as well as in frequency domain. All
post-processing related to the ACA solver is possible to be run directly after
step 2.

In the following these three steps are given in details.<br>
A detailed overview of all possible properties is given
[here](@ref propertyPageACA).

## Generating input # {#ugACA_generatingInput}

### Using sampled data from a previous solver run
To generate input data to be solved with the FWH method at first a FWH permeable
surface needs to be defined. This is done by providing a triangulated surface
geometry written in an ASCII .stl file format. These file may be generated using
a python script, paraview, or an external meshing tool.
It is possible that the surface is distributed among multiple files.
Assuming in the following these files are named `fwh0.stl`, `fwh1.stl` and so
forth.

To sample flow or acoustic field data in the previous solver run the
[post processing](@ref ugPostprocessing) needs to be enabled to activate the
surface sampling options
```PYTHON
  #..stuff relevant for the FV, LB, or DG solver run..

  # The following block is described in post processing page
  solvertype.1                = "MAIA_POST_DATA"
  geometryInputFileName.1     = "geometry.toml"
  outputDir.1                 = "out_surfaceSampling/"
  solvertype.1                = "POST_DATA"
  postprocessingOrder_0       = [0,1]
  noPostProcessing            = 1
  postProcessingSolverIds     = [0]
  postProcessingType_0        = "POSTPROCESSING_LB"

  # Enable the surface sampling post processing operation
  postProcessingOps_0         = ["PP_SURFACE_SAMPLING_IN"]

  # Configure sampling properties
  surfaceDataStartTimeStep    = 396800    # Start time step of sampling
  surfaceDataEndTimeStep      = 1984000   # End time step of sampling
  surfaceDataSampleInterval   = 64        # Interval of sampling data
  surfaceDataWriteInterval    = 9920      # Interval of writing data in a sampled file

  # Name of the surface/s to be considered for sampling
  surfaceDataFileName_0       = "stl/fwh0.stl"
  surfaceDataFileName_1       = "stl/fwh1.stl"
```
Here, every `9920` time steps for each of the two surfaces a surface sample file
is written. Each of these files contains
`surfaceDataWriteInterval/surfaceDataSampleInterval = 155` temporal samples for
each triangle defined in the `surfaceDataFileName` geometry. Thus, at the end
per each surface data
`(surfaceDataEndTimeStep-surfaceDataStartTimeStep)/surfaceDataWriteInterval=160`
sampled files are located in the directory `out_surfaceSampling/`.
The first file will be named with the suffix id number
`surfaceDataStartTimeStep/surfaceDataWriteInterval = 40`, i.e., the number
starts from the first time step even if later sampling start time step is chosen.

### Using analytical input function
The generation of analytical input data is done right in advance of a simulation
run. This is documented [here](@ref ugACA_setInputData_fromAnalytical) in the following
section.

## Far-field prediction # {#ugACA_farFieldPrediction}
To predict the acoustic far-field using the ACA solver three steps are relevant
for setting up a case
1. Choose the [numerical method](@ref ugACA_setNumericalMethods) see
[here](@ref nmACA) for more background information
2. Specify your [input data](@ref ugACA_setInputData) of the FWH permeable
surface
3. Set [boundary conditions](@ref ugACA_setBC), if necessary
4. Define the [observer points](@ref ugACA_setObserverPoints) of interest
5. (Choose relevant [post-processing](@ref ugACA_postProcessing) routine,
described in the next section)

In the following the different parts of a ACA property file are described
individually, whereas a merged file is given here in advance:
<details>
<summary>(Click here to expand complete property file)</summary>
More options and details about different types, e.g. for `windowType`, are
documented [here](@ref propertyPageACA).
```python
#-------------------------------------------------------------------------------
#---GENERAL---------------------------------------------------------------------
#-------------------------------------------------------------------------------
# This part deals with general properties described previously in the documention.

#---APPLICATION SETTINGS--------------------------------------------------------
flowSolver                        = True
gridGenerator                     = False

#---FLOW VARIABLES--------------------------------------------------------------
Ma                                = [0.1, 0.0, 0.0]

#---I/O-------------------------------------------------------------------------
outputDir.default                 = "./out_fwh/"

#---NUMERICAL PROPERTIES--------------------------------------------------------
nDim                              = 3
noSolvers                         = 1

executionRecipe                   = "RECIPE_BASE"
solvertype.default                = "MAIA_UNIFIED"

#---SOLUTION PROPERTIES---------------------------------------------------------
timeSteps                         = 1
solutionInterval                  = -1
residualInterval                  = -1

#---MEMORY----------------------------------------------------------------------
scratchSize                       = 25.0
maxNoCells                        = 1000000                 # required by m-AIA for total scratch size

#-------------------------------------------------------------------------------
#---ACA (FWH)-------------------------------------------------------------------
#-------------------------------------------------------------------------------
# This part deals with the properties relevant for the ACA solver using the FWH method

#---1. NUMERICAL PROPERTIES-----------------------------------------------------
solvertype.0                      = "MAIA_ACOUSTIC_ANALOGY" # Choosing the ACA solver
acousticMethod                    = "FWH"                   # Choosing the FWH method
transformationType                = 1                       # Here, 1: Fast-Fourier transformation
windowType                        = 3                       # Type of the window function. Here, 3: Hamming

#---2. Input data---------------------------------------------------------------
# To deal with input from different solvers a unit conversion step might be
# performed in advance of solving the FWH. The format of this conversion is
# defined by the acaResultsNondimMode property. Also, this property handles the
# form of the output units. For more reference it is referred to the property
# page of the ACA solver.
acaResultsNondimMode              = 2                       # Here, 2: LB input and Infty State output 

#--input: from sampling data
generateSurfaceData               = False                   # Disabling the analytical source term generation

useMergedInputFile                = False                   # If all sample files are merged together
inputDir                          = "./out_surfaceSampling/" # Directory where sample files are located
noSurfaces                        = 2                       # Number of FWH permeable surfaces
surfaceIds                        = [0,1]                   # Ids of the FWH permeable surface sample files

# NOTE:
# The number of samples should be a potency of 2. Otherwise no FFT can be used
# and computation is slowed down heavily.
# Here, each file contains 155 temporal samples
inputFileIndexStart               = [94]                    # Sample file index start
inputFileIndexEnd                 = [200]                   # Sample file index end
acaNoSamples                      = 16384                   # Number of temporal samples that shall be used

## Alternativ: Using a merged input file
#useMergedInputFile                = True                    # If all sample files are merged together
#inputDir                          = ""                      # Directory where sample files are located
#surfaceDataFileName_0             = "surfaceData.Netcdf"    # Name of the merged sample file

# #--input: from analytical function
# # Alternativ: Input is prescribed by analytical source term as discussed in
# # the previous section
# generateSurfaceData           = True                        # Enabling the analytical source term generation
# generateSurfaceDataNoSamples  = 64                          # Number of surface data samples to be generated
# #- Description of the analytical source term
# sourceType                    = 2                           # Type of source term (Here, 2 = Quadrupole)
# omega_factor                  = 2                           # Omega = pi * c_infty * omega_factor
# noPeriods                     = 8.3                         # Number of periods. Default = 1
# dimless                       = True                        # Generation of sound source without dimensions
# #- Description of the FWH permeable surfaces
# noSurfaces                    = 2                           # Number of FWH permeable surfaces to created
# # First FWH permeable surface
# surfaceType_0                 = 0                           # Geometry type of the surface (Here, 0 = plane)
# surfaceNormalDir_0            = 0                           # Surface normal dir (Here, 0 = -x direction)
# surfaceCoordinates_0          = [-3.0, -3.0, -3.0, -3.0, 3.0, 3.0 ] # Box of the plane [-x.-y.-z, +x,+y,+z]
# noSegments_0                  = [0, 40, 40]                 # Number segments in each Cartesian direction
# # Second FWH permeable surface
# surfaceType_1                 = 0
# surfaceNormalDir_1            = 1                           # Surface normal dir (Here, 1 = +x direction)
# surfaceCoordinates_1          = [3.0, -3.0, -3.0, 3.0, 3.0, 3.0]
# noSegments_1                  = [0, 40, 40]

#---3. Observer points----------------------------------------------------------
#--observer location: from .csv file
generateObservers                 = False                   # Generate observer (True) or use observerFile (False)
observerFileName                  = "observerLocation.csv"  # .csv file with observer point coordinates

# #--observer location: distribute observer points along a circular path
# # Alternative: Generate observer points, instead of using a .csv input file
# generateObservers                 = True
# observerRadius                    = 10.0                    # Radius of the circlular path
# noObservers                       = 100                     # Number of observer points generated

#---4. boundary condition---------------------------------------------------
# # Here, two symmetry planes are set, with normal in +y and in +z direction.
# symmetryOrigin = [0.0, 0.0, 0.0,                            # Origin of plane 1
#                   0.0, 0.0, 0.0 ]                           # Origin of plane 2
# symmetryNormal = [0.0, 0.0, 1.0,                            # Normal of plane 1
#                   0.0, 1.0, 0.0 ]                           # Normal of plane 2

#---5. ACA SPECIFIC POST PROCESSING---------------------------------------------
# The following is presented in the next section in more details
acaPostprocessingMode             = False                   # Apply after run (False), or on data already generated (True)
noPostprocessingOps               = 2                       # Number of post processings operations to be used
postprocessingOps                 = [0, 2]                  # Here, 0:RMS of pressure, 2: SPL
```
</details>

### Numerical methods # {#ugACA_setNumericalMethods}
```python
#---1. NUMERICAL PROPERTIES-----------------------------------------------------
solvertype.0                      = "MAIA_ACOUSTIC_ANALOGY" # Choosing the ACA solver
acousticMethod                    = "FWH"                   # Choosing the FWH method
transformationType                = 1                       # Here, 1: Fast-Fourier transformation
windowType                        = 3                       # Type of the window function. Here, 3: Hamming
```
The name of the solver is `MAIA_ACOUSTIC_ANALOGY`, which features an
`acousticMethod` in the moment named `"FWH"`
(see [mathematical model](@ref mmFWH) and [numerical method](@ref nmACA)).
As the FWH equation is solved in frequency domain and sampled input data and
some of the output data is in time domain a transformation (here Fourier
transformation) between these domains is required and set via
`transformationType`.
A Fourier transformation requires an infinite long temporal sample. Hence,
spectral leakage may occur when using a restricted input signal length. This is
reduced by choosing an appropriate windowing function via `windowType`.

### Input data of FWH permeable surface # {#ugACA_setInputData}
Data can be made non-dimensional in different ways. Therefore, the flag
`acaResultsNondimMode` is used to specify different modes of assuming input data
dimensions and writing output data with a certain dimension.
E.g., the mode `2` assumes that the input data is made non-dimensional as done in
the [LB solver](@ref LB_nonDim) while the output is made non-dimensional by the
infinity state as described in the [numerical method of the ACA solver](@ref nmACA_nonDim).
```python
#---2. Input data---------------------------------------------------------------
# To deal with input from different solvers a unit conversion step might be
# performed in advance of solving the FWH. The format of this conversion is
# defined by the acaResultsNondimMode property. Also, this property handles the
# form of the output units. For more reference it is referred to the property
# page of the ACA solver.
acaResultsNondimMode              = 2                       # Here, 2: LB input and Infty State output 
```
As already outlined in the [generating input](@ref ugACA_generatingInput)
section different ways of defining/reading input data are possible.

@anchor ugACA_setInputData_fromInputFile
To read **input from sampling files** obtained from a previous solver run as
described in [generating input](@ref ugACA_generatingInput) the following set of
properties is relevant
```python
#--input: from sampling data
generateSurfaceData               = False                   # Disabling the analytical source term generation

useMergedInputFile                = False                   # If all sample files are merged together
inputDir                          = "./out_surfaceSampling/" # Directory where sample files are located
noSurfaces                        = 2                       # Number of FWH permeable surfaces
surfaceIds                        = [0,1]                   # Ids of the FWH permeable surface sample files

# NOTE:
# The number of samples should be a potency of 2. Otherwise no FFT can be used
# and computation is slowed down heavily.
# Here, each file contains 155 temporal samples
inputFileIndexStart               = [94]                    # Sample file index start
inputFileIndexEnd                 = [200]                   # Sample file index end
acaNoSamples                      = 16384                   # Number of temporal samples that shall be used
```
In the case where multiple data on multiple surfaces is sampled, these are
selected using `noSurfaces` and the corresponding Id via `surfaceIds`.<br>
It is possible to use only a subset of the sampling files controlled via the
properties `inputFileIndexStart` and `inputFileIndexEnd`. Furthermore, the
number of temporal samples read from these files is controlled via
`acaNoSamples`.

@note
```python
useMergedInputFile = True
```
Using this property requires the usage of a python script located in the m-AIA
tools repository. With this it's possible to merge multiple sample files into a
single file in advance of an ACA solver run.

@anchor ugACA_setInputData_fromAnalytical
Input **from an analytical function** is possible by considering two steps.
First, a function describing the acoustic sources is to be selected.
Second, the shape of the used FWH permeable surface/s needs to be described.<br>
The source term type is choosen via `sourceType`, whereas it's frequency, signal
length, and units are controlled via the following properties
```python
#--input: from analytical function
# Alternativ: Input is prescribed by analytical source term as discussed in
# the previous section
generateSurfaceData           = True                        # Enabling the analytical source term generation
generateSurfaceDataNoSamples  = 64                          # Number of surface data samples to be generated
#- Description of the analytical source term
sourceType                    = 2                           # Type of source term (Here, 2 = Quadrupole)
omega_factor                  = 2                           # Omega = pi * c_infty * omega_factor
noPeriods                     = 8.3                         # Number of periods. Default = 1
dimless                       = True                        # Generation of sound source without dimensions
```
Multiple number of discrete FWH permeable surface/s might be defined using the
following
```python
#- Description of the FWH permeable surfaces
noSurfaces                    = 2                           # Number of FWH permeable surfaces to created
# First FWH permeable surface
surfaceType_0                 = 0                           # Geometry type of the surface (Here, 0 = plane)
surfaceNormalDir_0            = 0                           # Surface normal dir (Here, 0 = -x direction)
surfaceCoordinates_0          = [-3.0, -3.0, -3.0, -3.0, 3.0, 3.0 ] # Box of the plane [-x.-y.-z, +x,+y,+z]
noSegments_0                  = [0, 40, 40]                 # Number segments in each Cartesian direction
# Second FWH permeable surface
surfaceType_1                 = 0
surfaceNormalDir_1            = 1                           # Surface normal dir (Here, 1 = +x direction)
surfaceCoordinates_1          = [3.0, -3.0, -3.0, 3.0, 3.0, 3.0]
noSegments_1                  = [0, 40, 40]
```

### Boundary conditions # {#ugACA_setBC}
For certain problems the usage of a symmetry boundary might be useful if the
sampled data is not available. Therefore, a symmetry plane is defined by a point
(`symmetryOrigin`) and a surface normal (`symmetryNormal`). For usage of
multiple planes these two properties are simply continued
```python
# Here, two symmetry planes are set, with normal in +y and in +z direction.
symmetryOrigin = [0.0, 0.0, 0.0,                            # Origin of plane 1
                  0.0, 0.0, 0.0 ]                           # Origin of plane 2
symmetryNormal = [0.0, 0.0, 1.0,                            # Normal of plane 1
                  0.0, 1.0, 0.0 ]                           # Normal of plane 2
```

### Definition of observer points # {#ugACA_setObserverPoints}
The location of the observer points is described either via a `.csv` file or by
using an analytical description, flagged by the property `generateObservers`.<br>
The location of the `.csv` file is given with
```python
#--observer location: from .csv file
generateObservers                 = False                   # Generate observer (True) or use observerFile (False)
observerFileName                  = "observerLocation.csv"  # .csv file with observer point coordinates
```
Here, `observerLocation.csv` might look as following, defining two arbitrarily
located observer points
```
1.00001 10.0  -2.0
20.0    10.0  -2.0
```
To generate observer points distributed along a circular path an
`observerRadius` as well as a `noObservers` needs to be given
```python
#--observer location: distribute observer points along a circular path
# Alternative: Generate observer points, instead of using a .csv input file
generateObservers                 = True
observerRadius                    = 10.0                    # Radius of the circlular path
noObservers                       = 100                     # Number of observer points generated
```

## Post-processing # {#ugACA_postProcessing}
Different post processing operations are provided by the ACA solver, which might
be applied directly after the solver run, or in a post-step on observer data
already generated using the post processing mode
(`acaPostprocessingMode = True`). The number of post processing operations
(`noPostprocessingOps`) and the type of the `postprocessingOps` is given here
for an example
```python
#---5. ACA SPECIFIC POST PROCESSING---------------------------------------------
# The following is presented in the next section in more details
acaPostprocessingMode             = False                   # Apply after run (False), or on data already generated (True)
noPostprocessingOps               = 2                       # Number of post processings operations to be used
postprocessingOps                 = [0, 2]                  # Here, 0:RMS of pressure, 2: SPL
```
all operations available are documented in the
[ACA property page](@ref propertyPageACA).
