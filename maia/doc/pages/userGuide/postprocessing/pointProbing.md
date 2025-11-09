# Point Probing # {#ppPointProbing}

Point probing allows you to specify a point coordinates and generate output files containing
the variables of the cells at the point coordinates (nearest neighbor).
The point coordinate is specified by
```
pp_probeCoordinates
```
The number of points is determined by the given total number of points.
The variables at the probe points are saved in ```probe_<probepointId>.dat``` in the folder
```
pp_probePath
```
The default folder is ```./probes/```

The operation is triggered by setting one of the following properties
```
postProcessingOps_<postprocessinginstanceId> = [PP_PROBE_POINT_PRE]
postProcessingOps_<postprocessinginstanceId> = [PP_PROBE_POINT_IN]
postProcessingOps_<postprocessinginstanceId> = [PP_PROBE_POINT_POST]
```
