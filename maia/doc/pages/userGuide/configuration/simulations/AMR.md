# Adaptive mesh refinement (AMR) # {#ugAMR}

## General

When using m-AIA with a Cartesian solver, you have the possibility to use solution adaptive mesh refinement. This allows you to refine or coarsen the grid during the simulation based on the solution state.
To active adaptive mesh refinement for your simulation you simply set

```
adaptation = true
```
To control the time step interval at which the adaptive mesh refinement is performed, you have to set

```
adaptationInterval = 42
```
Here, the adaptation is performed every 42 time steps, excluding the initial time step. Sometimes if you don't want to adapt your grid based on initial fluctuations in your solution, it is a good idea to delay the start of the AMR.
This can be done by setting

```
adaptationStart = 200
```
If you want to adapt you grid based on the initial condition of your simulation, you can do so by using

```
initialAdpatation = true
```
In general, the AMR is happening between the maximum uniform refinement (```maxUniformRefinementLevel```) level and the maximum refinement level (```maxRfnmtLvl```). In other word, no cell can be coarsened below the uniform grid or refined above the maximum refinement level.
For specific applications, it can be of interest to increase the maximum uniform or maximum refinement level, 
which is also possibly.

Additionally to regular adaptation intervals, AMR can be triggered based on the solution evolution, 
i.e. when embedded bodies or Lagrangian particles move in the proximity of their refinement.
This feature is activated by default when using the level-set or Lagrangian particle solver.

An extensive list of all properties related to AMR can be found [here.](@ref propertiesAMR)

## Sensors

During the adaptation the algorithm has to decide whether to **refine**, **coarsen** or **keep** the grid cell as it is. This is done by using solution based sensors which represent different metrics associated with the local resolution of the grid.
For example, if your solution method should be able to capture steep gradients, a certain local grid resolution is needed. Using the gradient of your solution as sensor metric, the AMR will then ensure a fine grind in regions of high gradients.
When using multiple sensors, which might yield contradictory results, **refinement** will always prevail.
There are two categories of sensors:

@paragraph sensorsIndirect Indirect sensors

These sensors a evaluated and return a floating point number for each cell of the grid (e.g. gradient, vorticity).
This metric is then evaluated statistically, such that cells with a sensor value sufficiently above/below the root-mean-square (RMS) will be refined/coarsened.
All remaining cells are kept at their current level. To specify the ratio between the refinement threshold and the RMS value of the sensor you have to set ```sensorWeight``` for each sensor.
The coarsening threshold is controlled via the global property ```coarseRatio``` which is the ratio between the refinement and coarsening threshold.

@paragraph sensorsDirect Direct sensors

These sensors enforce the adaptation of the grid directly without comparing the cell values to each other. For example, if there is a cell type in the simulation which always needs to be refined, no comparison is needed. For this type of sensor the weight has to be set to -1.

@note Which sensors are available for your simulation depends on the specific solver used. Please refer to the user guide of the solver for a list of available adaptation sensors.

### Example

In this example we want to adaptively refine the grid during the simulation of a moving body in a fluid. We want to make sure, that the interface between our fluid and solid phase is always refined. Additionally we want to resolve the body wake using the divergence and vorticity as indicator. The adaptation should be performed every 100 time steps, starting with the initial condition. The corresponding properties would look like this:

```
adaptation = true
adaptationStart = 0
adaptationInterval = 100
initialAdaptation = true

sensorType = ["INTERFACE", "DIVERGENCE", "VORTICITY"]
sensorWeight = [-1, 1.5, 1.5]
coarseRatio = 0.2
```

### Advanced features

If you want to visualize the sensor values directly, you can set ```saveSensorData = true```. This will write a separate output file each time the adaptation is triggered.

To restrict the AMR to refining cells only, i.e., never coarsen cells even if the sensor value is below the threshold, you could set ```allowCoarsening = false```.
