# Postprocessing # {#ugPostprocessing}
[TOC]

## Overview
When using m-AIA with a Cartesian solver, you have the possibility to use Postprocessing.

The postprocessing is activated by setting:
```
postProcessing = true
```
Multiple postprocessing instances can be created to do Post-Processing for a multi-solver simulation.
The number of postprocessing instances is specifies by the property
```
noPostProcessing
```
The assignment of a postprocessing instance with the corresponding physics solver is done by the property
```
postProcessingSolverIds
```
The postprocessing type specifying the type of physics solver and the postprocessing operations of a postprocessing instance is specified by the property
```
postProcessingType
```
Possible Post-Processing types are
```
"POSTPROCESSING_FV"
"POSTPROCESSING_LS"
"POSTPROCESSING_LB"  
"POSTPROCESSING_DG"
"POSTPROCESSING_FVLPT"  
"POSTPROCESSING_LBLPT"
```   

The Post-Processing operations are specified by
```
postProcessingOps
```
The list of Post-Processing operations is given below.

@subpage ppAverage  
@subpage ppMovingAverage    
@subpage ppPointProbing  
@subpage ppLineProbing   
@subpage ppArbLineProbing   
@subpage ppSliceProbing   
@subpage ppArbSliceProbing  
@subpage ppRecuceLevel  
@subpage ppSpatialAverage

Every Postprocessing functionality can be labeled in three different ways   
`PRE` PreSolve: Call of postprocessing operations before simulation is started (performs postprocessing on restartFiles that are loaded)   
`IN` InSolve: Call of postprocessing operations during simulation   
`POST` PostSolve: Call of postprocessing operations after simulation (i.e. probing on averaged variables)

Here the number of the postprocessing instance is specified by appending ```_<postprocessinginstanceId>```.

Depending on whether an additional postdata solver is called, the Post-Processing practice can be classified into two general types:
<ol>
<li>Post-Processing with an additional Postdata solver</li>
<li>Post-Processing without the Postdata solver for special IO</li>
</ol>

##Post-Processing with an additional Postdata solver
The first type of Post-Processing calculates additional flow properties
(i.e. temporal average, moving average...), which are stored in a separate output file.
For that a seperate solver has to be created in the grid file, see [Mesh Generation](@ref ugMeshGeneration).
This allows the Postdata solver to be a subset of the global grid, such that postprocessing operations can be
performed on a subset of the global grid (i.e. calculation of average in a specific region).
Ideally, the postdata solver should be the solver with the highest solverId.

Additionally, the following properties have to be specified to configure the Postdata solver:   
```
solvertype.<postdatasolverId> = "MAIA_POST_DATA"
noVariables.<postdatasolverId> = 5   # for example 
```
The variable ```postdatasolverId``` is set to the solver-ID of the Postdata solver in the grid file.

Make sure, that ```multiSolverGrid = true```.

##Post-Processing for special IO
The second type allows the output of reduced solver data (i.e. line probing, slice probing...).
For that, a separate Postdata solver is not required, unless the probing is performed on averaged variables (PostSolve).

##Working Example
The following example shows how to set up the Post-Processing instance for a coupled FV-LB simulation (FV-Solver solverId = 0, LB-Solver solverId = 1, Postdata solverId = 2).
```
 postProcessing          = true
 noPostProcessing        = 2
 postProcessingSolverIds = [0,1]
 postProcessingType_0    = ["POSTPROCESSING_FV"]
 postProcessingOps_0     = ["PP_AVERAGE_IN", "PP_PROBE_LINE_IN"]
 postProcessingType_1    = ["POSTPROCESSING_LB"]
 postProcessingOps_1     = ["PP_AVERAGE_IN"]
 
 solvertype.2            = "MAIA_POST_DATA"
 # number of variables (5 primitive variables + 7 additional variables for pp_square calculation)
 noVariables.2           = 12 
 
 # Properties for PP_AVERAGE_IN
 pp_square               = true
 pp_averageStartTimestep = 1
 pp_averageStopTimestep  = 10
 pp_averageInterval      = 1
 pp_averageRestart       = false

 # Properties for PP_PROBE_LINE_IN
 pp_probeLineDirection   = [0, 2]
 pp_probeLineCoordinates = [-0.1, -0.1, 0.0, 0.0]
```
