# Spatial Average # {#ppSpatialAverage}

This Post-Processing operation performs a spatial average in Cartesian direction. Set
```
pp_spatialDirection1
pp_spatialDirection2
```
If `pp_spatialDirection1` and `pp_spatialDirection2` are both -1, spatial averaging of the
whole flow field on a single point is performed. If one of the two directions is unequal -1,
spatial averaging on a line in the specified direction is performed. If both directions are
unequal -1 and not equal to each other, same is done for a slice in these two directions (3D only).

The operation is triggered by setting one of the following properties
```
postProcessingOps_<postprocessinginstanceId> = [PP_SPATIAL_AVERAGE_PRE]
postProcessingOps_<postprocessinginstanceId> = [PP_SPATIAL_AVERAGE_IN]
postProcessingOps_<postprocessinginstanceId> = [PP_SPATIAL_AVERAGE_POST]
```
The ```PRE``` and ```POST``` modes require restart files of the averaging timesteps.   
The ```IN``` mode performs averaging done during the solver run, which does not require restartfiles.   
For that, set
```
MInt pp_averageStartTimestep
MInt pp_averageStopTimestep
MInt pp_averageInterval
```
