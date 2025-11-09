# Line Probing # {#ppLineProbing}

Line probing allows you to specify multiple lines pointing in a cartesian coordinate direction
and generate output files containing the variables of the cells intersected by each line
(nearest neighbor). The cartesian direction is set by the property ```pp_probeLineDirection```, hereby is
```
pp_probeLineDirection = 0 -> x-direction
pp_probeLineDirection = 1 -> y-direction
pp_probeLineDirection = 2 -> z-direction
```
The position of the line is defined by ```nDim-1``` points set in ```pp_probeLineCoordinates```,
which define one point on the line in the remaining directions.
The number of lines is determined by the given total number of points.

The operation is triggered by setting one of the following properties
```
postProcessingOps_<postprocessinginstanceId> = [PP_PROBE_ARB_LINE_PRE]
postProcessingOps_<postprocessinginstanceId> = [PP_PROBE_ARB_LINE_IN]
postProcessingOps_<postprocessinginstanceId> = [PP_PROBE_ARB_LINE_POST]
```

Additionally, the following properties have to be set:
```
pp_probeLineInterval
pp_probeLineStartTimestep
pp_probeLineStopTimestep
```

##Working example
The following example loads a mean file and performs a line probing on the averaged variables
for two lines. The first line runs along the x-direction and is located at (y=0.1, z=0.1).
```
### Postprocessing instance
postProcessing          = true
noPostProcessing        = 1
postProcessingSolverIds = [0]
postProcessingType_0    = "POSTPROCESSING_FV"
postProcessingOps_0     = ["PP_PROBE_LINE_POST"]

### Postdata
solvertype.1  = "MAIA_POST_DATA"
noVariables.1 = 5

### Properties for PP_PROBE_LINE_POST
pp_fileName             = "out/Mean_s1_00002005-00002050.Netcdf"
pp_averageStartTimestep = 2005
pp_averageStopTimestep  = 2050
pp_averageInterval      = 5
pp_probeLineDirection   = [0]
pp_probeLineCoordinates = [0.1, 0.1]
```
