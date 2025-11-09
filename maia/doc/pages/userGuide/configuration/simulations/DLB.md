# Dynamic Load Balancing (DLB) # {#ugDLB}

## General

To achieve a high parallel efficiency, m-AIA incooperates a dynamic load balancing approach.
When a dynamic load balancing (DLB) is called, the cartesian cells 
and their associated information are re-distributed among all compute cores.
This is especially important for applications with a time varying computational load, 
i.e. due to adaptive mesh refinement or differing computing architecture.
To activate the dynamic load balancing simply set
```
balance = true
```

For dynamic load balancing on static grids, the following properties
```
loadBalancingInterval = 200
loadBalancingStopTimeStep = 2000
```
are relevent to trigger the interval and a stopping time step for the DLB.

On a cartesian grid with AMR, the computational load may vary strongly after an adaptation time step.
Thus, it is meaningful to trigger DLB not based on time steps, but based on adaptation itself.
This can be controlled by the following properties:
```
balanceAfterAdaptation = true
balanceAfterAdaptationInterval = 2
```
which ensure a balance after every second adaptation call.

### Cell Weight Computation
A main criterion for the DLB is the type of cell weight computation.
In m-AIA, the cell weight can be either defined as a static cell weight 
or be computed dynamically based on domain run- and idle timers.
This distinction is controlled by the ```loadBalancingMode``` property, which can be set to:  
0: for a static weight computation  
1: for a dynamic weight computation.

The dynamic weight computation requires little user input and computes cell weights for 
different solver cell types based on the run timers of all domains 
and the distribution of the different cell types among domains.

For static weight computation values for the different cell weights must be defined a priori 
by the user. Depending on the application and multi-physics problem, 
different solver- and cell specific weights are optimal. 
However, default values are implemented and can be tuned 
for the specific application by setting the properties listed here.

### Partitioning

The distribution of the computational load among all compute cores, 
is tied strongly to the partitioning algorithm. 
The following properties can be used to control the partition-cell generation and distribution.
```
dlbUpdatePartitionCells = true
partitionCellOffspringThreshold = 5000
partitionCellWorkloadThreshold  = 25000
dlbPartitionMethod    = "DLB_PARTITION_WEIGHT"
```

An extensive list of all DLB related properties can be found [here.](@propertiesDLB)
