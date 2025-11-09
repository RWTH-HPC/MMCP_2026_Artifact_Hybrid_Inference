# Abritrary Slice Probing # {#ppArbSliceProbing}

The abritrary slice probing outputs the postprocessing variables on an abritrary slice that
is defined by three points specified by the property ```pp_arbSlicePoints```.
The number of slices is determined by the given total
number of points. Given points a, b and c for a slice (with vectors b-a and c-a nonparallel),
the positions on the slice are determined by p_ij = a + 1/(i-1)*(b-a) + 1/(j-1)*(c-a).
The parameters i and j (number of points in direction b-a and c-a) are determined by the
property ```pp_noPointsArbSlice```.  
The output is a single output file named ```probeArbitrarySlices_<globalTimeStep>.Netcdf``` containing all slices.

The operation is triggered by setting one of the following properties
```
postProcessingOps_<postprocessinginstanceId> = [PP_PROBE_ARB_SLICE_PRE]
postProcessingOps_<postprocessinginstanceId> = [PP_PROBE_ARB_SLICE_IN]
postProcessingOps_<postprocessinginstanceId> = [PP_PROBE_ARB_SLICE_POST]
```

Additionally, the properties
```
MInt pp_averageStartTimestep
MInt pp_averageStopTimestep
MInt pp_averageInterval
```
have to be set.
