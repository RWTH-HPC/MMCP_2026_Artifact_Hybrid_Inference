# Moving Average # {#ppMovingAverage}

This Post-Processing operation performs a moving average on a rectangular time window.
The time window is set by the property by
```
pp_movingAverageDataPoints
```

The moving average starts at
```
pp_movingAverageStartTimestep
```
and ends at
```
pp_movingAverageStopTimestep
```
with an interval
```
pp_movingAverageInterval
```

The operation is triggered by setting one of the following properties
```
postProcessingOps_<postprocessinginstanceId> = [PP_MOVING_AVERAGE_PRE]
postProcessingOps_<postprocessinginstanceId> = [PP_MOVING_AVERAGE_IN]
postProcessingOps_<postprocessinginstanceId> = [PP_MOVING_AVERAGE_POST]
```
