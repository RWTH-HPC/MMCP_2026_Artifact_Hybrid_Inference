# Arbitrary Line Probing # {#ppArbLineProbing}

The arbitrary line probing outputs the postprocessing variables on an arbitrary line that is defined by two points specified by the properties ```pp_arbLinePoints```, which defines the start and end point of the line. The number of lines is determined by the given total number of points. The arbitrary line is divided into equidistant points triggered by the property ```pp_noArbLineIds```. The postprocessing variables are interpolated onto the equidistant points by the least square method. The output is saved as a Netcdf-file with the name ```probeArbitraryLines_<globalTimeStep>.Netcdf```.

The operation is triggered by setting one of the following properties
```
postProcessingOps_<postprocessinginstanceId> = [PP_PROBE_LINE_PRE]
postProcessingOps_<postprocessinginstanceId> = [PP_PROBE_LINE_IN]
postProcessingOps_<postprocessinginstanceId> = [PP_PROBE_LINE_POST]
```

Additionally, the properties
```
MInt pp_averageStartTimestep
MInt pp_averageStopTimestep
MInt pp_averageInterval
```
have to be set.
