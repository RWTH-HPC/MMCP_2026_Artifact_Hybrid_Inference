# Parallelization # {#ugParallelization}

## Partitioning

In principal, the partitioning of the grid occures on the `minLevel` of the cartesian grid. 
However, for applications with large and local refinement run on many compute-cores, 
it is meaningful to create additional partition cells on higher refinement levels. 
This process is called "partition level shift" and the parents of the partition cells are called "partition level anchestors". 
The introduction of these additional partition cells is trigged by the following two properties:
```
partitionCellOffspringThreshold  = 512
partitionCellWorkloadThreshold = 1024
```
This indicates, that a minLevel cell with more than 512 childs 
or with a combined workload of all children of more than 1024 becomes a new partition cell.

Note, that the additional partition cells allow for a better load distribution for high-performance computing, but also increase the comunication between different cores.

In general, the partitioning of the cartesian grid is performed at the beginning of the simulation 
and a re-partitioning is performed at every restart, unless 
```
updatePartitionCellsOnRestart = false
``` 

Updating a partition is especially important for dynamically changing grids, 
i.e. grids with adaptive mesh refinement.

In general m-AIA allows for a an MPI and hybrid MPI and open-MP parallelization. 
For the MPI comminucation of the cartesian grid, additioan comminucation cells are created at the boundary of two MPI ranks. These cells are calles `halo-Cells`. Halo-Cells are required for the 
reconstruction of slopes or fluxes between the different MPI ranks. 
The number halo layers can be triggered with the property `noHaloLayers` and is dependend on the solver solution scheme. In gerneral 2 halo-layers are sufficient for FV-applications.


## Memory

For complex large-scale simulations with multiple coupled solvers it is
meaningful to have a look at the total allocated memory, e.g., to avoid
out-of-memory crashes or unnecessary excessive memory allocations which might
impair the overall performance of the simulation.

The size of the pre-allocated memory can be controlled mainly via:
```
maxNoCells = <N>
scratchSize = <X>
```
with `maxNoCells` the maximum number of cells of the grid on a
single rank and `scratchSize` determining the total size of the
`scratch`-buffers used for large temporary memory allocations.
The total scratch memory allocation is `sizeof(MFloat) * scratchSize * maxNoCells`.

The memory statistics output can be controlled via the property:
```
displayMemoryStatistics = true
```

Exemplary memory statistics output of a benchmark on 4096 nodes of the HAWK
system:
```
******************************* MEMORY STATISTICS *******************************
***** Comment: After run loop - #ranks: 524288
***** Location: void MAIAApplication::run() [with int nDim = 3] (/zhome/academic/HLRS/xac/xacanie/code_maia/MAIA-master/src/application.cpp:1131)
*****                                                                               // Memory statistics: physical=RAM vs. total allocation size
***** Average memory usage: physical = 1390.31 MB; allocation = 1646.25 MB          // Average memory per rank
***** Minimun memory usage: physical = 1370.01 MB; allocation = 1637.69 MB          // Minimum memory over all ranks
***** Maximum memory usage: physical = 3199.09 MB; allocation = 3376.36 MB          // Maximum memory over all ranks
***** Maximum diff in memory usage: physical = 2.97266 MB; allocation = 1.93359 MB  // Maximum difference in memory usage compared to last `memory statistics` evaluation
***** Total physical memory usage (RAM): 711837 GB                                  // Total RAM usage
***** Diff total physical memory usage (RAM): 423.216 GB                            // Difference in RAM usage compared to last evaluation
***** Total allocation size (Virtual Memory): 842883 GB                             // Total allocation size
***** Diff total allocation size (Virtual Memory): 68.5437 GB                       // Difference in total allocation size compared to last evaluation
*****
***** Maximum stack memory: 8220 KB; stack limit 1.80144e+16 KB                     // Maximum and limit of stack memory usage
*****
***** Minimum available memory per node (meminfo): 86.8946 GB                       // Minimum still available memory per node
***** Minimum free memory per node (RAM): 87.5009 GB
******************************* MEMORY STATISTICS *******************************
```


## Parallel efficiency


### Performance evaluation

There are different ways to evaluate the overall performance of a simulation.

Setting the property
```
aliveInterval = <N>
```
will enable a repeated output of the measured total runtime for each `N` simulation time steps.
This is a quick way to check/estimate how much time is required to perform a simulation and to detect odd performance behavious over the course of a simulation.

If load balancing is not enabled the DLB timers can be used for a performance output in an automatically determined interval based on the total number of time steps:
```
performanceOutput = <true/false> // enabled by default
```

Examplary performance output of a coupled FV-DG simulation run using 16384 MPI processes:
```
<m d="0" >456400 * Average time per step: 0.153592 ; local: 0.153427, idle/comp = 0.126337, idle/timePerStep = 0.112045
</m>
<m d="0" >456400 * Relative idle time: max = 0.625915, min 0.005925
</m>
<m d="0" >456400 * maxRunTime 0.152596
</m>
<m d="0" >456400 * Imbalance percentage: 25.56%, t_avg = 0.1136, t_max = 0.152608
456400 * Loads: max = 1.34338, min = 0.502796

 |------------------------------------------------------------------------------------------------------------------------------|
 | Load distribution at global timestep 456400     - imbalance 25.56% - t_avg 1.136e-01 - t_max 1.526e-01 - loads [0.503,1.343] |
 |------------------------------------------------------------------------------------------------------------------------------|
 | Load distribution at global timestep 456200     - imbalance 28.85% - t_avg 1.117e-01 - t_max 1.570e-01 - loads [0.511,1.405] |
 |------------------|------------------------------------------|------------------------------------------|---------------------|
 |     load bin     | current distribution - timestep   456400 | previous distribution - timestep  456200 |   curr/prev count   |
 |------------------|------------------------------------------|------------------------------------------|---------------------|
 | [ 1.450, 1.500 ) |                                          |                                          |        0 |        0 |
 | [ 1.400, 1.450 ) |                                          | @                                        |        0 |        2 |
 | [ 1.350, 1.400 ) |                                          | @                                        |        0 |        3 |
 | [ 1.300, 1.350 ) | #                                        | @                                        |       28 |       18 |
 | [ 1.250, 1.300 ) | #                                        | @                                        |       38 |       45 |
 | [ 1.200, 1.250 ) | #                                        | @                                        |       91 |       97 |
 | [ 1.150, 1.200 ) | ##                                       | @@                                       |      247 |      264 |
 | [ 1.100, 1.150 ) | #######                                  | @@@@@@@@@@                               |      878 |     1282 |
 | [ 1.050, 1.100 ) | ###########################              | @@@@@@@@@@@@@@@@@@@@@@@@@@@@@            |     3266 |     3503 |
 | [ 1.000, 1.050 ) |_########################################_|_@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@       _|     4831 |     4043 |
 | [ 0.950, 1.000 ) | ###########################              | @@@@@@@@@@@@@@@@@@@@@@@@                 |     3376 |     2983 |
 | [ 0.900, 0.950 ) | #############                            | @@@@@@@@@@@@@@@                          |     1639 |     1823 |
 | [ 0.850, 0.900 ) | #######                                  | @@@@@@@@@                                |      918 |     1110 |
 | [ 0.800, 0.850 ) | ####                                     | @@@@@                                    |      541 |      647 |
 | [ 0.750, 0.800 ) | ##                                       | @@                                       |      297 |      343 |
 | [ 0.700, 0.750 ) | #                                        | @                                        |      131 |      113 |
 | [ 0.650, 0.700 ) | #                                        | @                                        |       59 |       62 |
 | [ 0.600, 0.650 ) | #                                        | @                                        |       37 |       33 |
 | [ 0.550, 0.600 ) | #                                        | @                                        |        5 |        6 |
 | [ 0.500, 0.550 ) | #                                        | @                                        |        2 |        7 |
 | [ 0.450, 0.500 ) |                                          |                                          |        0 |        0 |
 |------------------|------------------------------------------|------------------------------------------|---------------------|
```


### Timings

Note: timings reported in `m_log` are specific to the global root domain (rank 0) and not meaningful for all ranks.

During a simulation timings of all solvers and couplers on all ranks can be collected and stored for later performance analysis:
```
writeSolverTimings = true            // Enable collection of timings
solverTimingsWriteInterval = <N>     // Write timings file every N time steps
solverTimingsSampleInterval = <M>    // Sample timings every M-th time step
writeAllSolverTimings = <true/false> // Output all timings (true) or a reduced set of the most essential timings (false)
```
The timings files can be postprocessed using a script in the Tools repository: postprocessing/process_block_timings.py


## TODO
* Change Initial Cell Weights (on static grids)
* partition parallel split
* Dynamic Load Balancing
* IO
* other best practices?
* auto-save feature?
