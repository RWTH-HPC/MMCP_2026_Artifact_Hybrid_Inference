# Slice Probing # {#ppSliceProbing}

Slice probing allows you to specify a plane spanning in two cartesian coordinate directions
and generate output files containing the variables of the cells intersected by the plane (nearest neighbor). In case the plane lies between cells, the cells in the slice do not have to be direct neighbors.
The cartesian direction is set by the property ```pp_probeSliceDir```, hereby is
```
pp_probeSliceDir = 0 -> x-direction
pp_probeSliceDir = 1 -> y-direction
pp_probeSliceDir = 2 -> z-direction
```
The position of the slice is determined by the position in slice direction with the property ```pp_probeSliceCoordinate```. The number of slices is determined by the given total number of points.

The operation is triggered by setting one of the following properties
```
postProcessingOps_<postprocessinginstanceId> = [PP_SLICE_LINE_PRE]
postProcessingOps_<postprocessinginstanceId> = [PP_SLICE_LINE_IN]
postProcessingOps_<postprocessinginstanceId> = [PP_SLICE_LINE_POST]
```

Additionally, the following properties have to be set:
```
pp_probeSliceInterval
pp_probeSliceStartTimestep
pp_probeSliceStopTimestep
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
postProcessingOps_0 = ["PP_PROBE_SLICE_POST"]

### Postdata
solvertype.1 = "MAIA_POST_DATA"
noVariables.1 = 18

### Properties for PP_PROBE_SLICE_POST
pp_fileName             = "out/Mean_s1_00002005-00002050.Netcdf"
pp_averageStartTimestep = 2005
pp_averageStopTimestep  = 2050
pp_averageInterval      = 5
pp_probeLineDirection   = [0]
pp_probeLineCoordinates = [0.1, 0.1]
```
