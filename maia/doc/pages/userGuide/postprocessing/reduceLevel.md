# Reduce Level # {#ppRecuceLevel}

This operation allows your to generate a coarser solution to you existing mesh.
The coarser level is set by
```
<MInt> pp_reductionLevel
```

The operation is triggered by setting one of the following properties
```
postProcessingOps_<postprocessinginstanceId> = [PP_REDUCE_TO_LEVEL_PRE]
postProcessingOps_<postprocessinginstanceId> = [PP_REDUCE_TO_LEVEL_POST]
```

You can also reduce the mesh resolution of a average solution. The average solution is
set by:
```
<MString> pp_ReStressesAverageFileName
```
This operation can be triggered by setting
```
postProcessingOps_<postprocessinginstanceId> = [PP_REDUCE_TO_LEVEL_AVERAGES_PRE]
```