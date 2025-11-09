# Average # {#ppAverage}

This Post-Processing operation performs a temporal average on the primitive variables of the assigned solver starting at the timestep ```pp_averageStartTimestep``` until ```pp_averageStopTimestep``` with an interval of ```pp_averageInterval```.   
The operation is triggered by setting one of the following properties
```
postProcessingOps_<postprocessinginstanceId> = [PP_AVERAGE_PRE]
postProcessingOps_<postprocessinginstanceId> = [PP_AVERAGE_IN]
postProcessingOps_<postprocessinginstanceId> = [PP_AVERAGE_POST]
```
The ```PRE``` and ```POST``` modes require restart files of the averaging timesteps.   
The ```IN``` mode performs averaging done during the solver run, which does not require restartfiles.

Without setting any other property, only the temporal average is calculated for the primitive variables.
The average variables will be stored in ```restartFile_<postdataId>_<globaltimestep>.Netcdf```, which will be named ```"um", "vm", "wm"(3D), "rhom", "pm"```.

##Additional calculations
```
<MBool> pp_square = true
```
This property triggers the calculation of the Reynold stresses   
@f$\left(\begin{array}{ccc}
\overline{u^\prime u^\prime} & \overline{u^\prime v^\prime} & \overline{u^\prime w^\prime} \\
\overline{v^\prime u^\prime} & \overline{v^\prime v^\prime} & \overline{v^\prime w^\prime} \\
\overline{w^\prime u^\prime} & \overline{w^\prime v^\prime} &\overline{w^\prime w^\prime}
\end{array}\right)@f$

as well as @f$\overline{p^\prime p^\prime}@f$.   
Increase the number of allocated variables of the postData (```noVariables.<postdataId>```) solver by 7.

```
<MBool> pp_skewness = true
```
This property determines if skewness is computed.   
Increase the number of allocated variables of the postData (```noVariables.<postdataId>```) solver by 3.

```
<MBool> pp_kurtosis = true
```
This property determines if kurtosis (and skewness) is computed.   
Increase the number of allocated variables of the postData (```noVariables.<postdataId>```) solver by 3.

```
<MBool> pp_averageVorticity = true
```
This property determines if the vorticity is included in the calculation.    
Increase the number of allocated variables of the postData (```noVariables.<postdataId>```) solver by ```nDim * 2 - 3```.

```
<MBool> pp_averageSpeedOfSound = true
```
Increase the number of allocated variables of the postData (```noVariables.<postdataId>```) solver by ```1 + nDim```.

##Special summation
In order to reduce rounding/truncation errors either activate kahan summation, which requires more memory,
or two-pass mode, with a higher time cost.
Two-pass is not applicable for the IN mode, since mean values are required for computation.
```
<MBool> pp_twopass = true
<MBool> pp_kahan = true
```

## Restarting Averging
Set the following properties as shown to restart the temporal average:
```
pp_averageRestart = true
restartFile.<postdataId> = true
```

##Working example
The following example restarts a temporal average calculation of a FV-Solver. In this example the ```restartTimeStep``` is between 100000 and 200000.
```
### Postprocessing instance
postProcessing          = true
noPostProcessing        = 1
postProcessingSolverIds = [0]
postProcessingType_0    = "POSTPROCESSING_FV"
postProcessingOps_0     = ["PP_AVERAGE_IN"]

### Postdata
solvertype.1  = "MAIA_POST_DATA"
# 5 primitive variables + 7 square variables + 3 vorticity variables
noVariables.1 = 15
restartFile.1 = true

### Properties for PP_AVERAGE_IN
pp_averageStartTimestep = 100000
pp_averageStopTimestep  = 200000
pp_averageInterval      = 5
pp_averageRestart       = true
pp_square               = true
pp_averageVorticity     = true
```
