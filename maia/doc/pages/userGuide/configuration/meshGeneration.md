# Mesh Generation # {#ugMeshGeneration} 
The m-AIA framework has its own tool for the parallel generation of hierarchical [Cartesian meshes](#nmCartesianGrid)
 which are used by the [FV](#nmFV), [FVMB](#nmFVMB), [LB](#nmLB), [DG](#nmDG), and the [LS](#mmLS) solver. Like the 
 m-AIA solver environment the grid generator receives the required information from a `property.toml`-file. To tell 
 m-AIA to run the grid generator instead of the solver environment set the following properties to:
```
##### properties.toml #####
gridGenerator = true  
flowSolver = false
```
In general, it is wise to use specific property files for the grid generator and the solver setup to avoid an 
excessive amount of properties that make it difficult to create a concise setup. In the following we will go 
through the various functionalities of the grid generator and the corresponding properties. 

## Setting input/output paths
To define the directory the generated grid is written to the property `outputDir` can be set accordingly. The name of 
the grid file can be defined setting the property `gridOutputFileName`. The name of the `geometry.toml`-file containing 
information on the files providing the geometry itself (and the according boundary conditions) is given by
`geometryInputFileName` Note that all file paths refer to the base directory of the setup. An examplary setup for these 
properties could be:
```
##### properties.toml #####
outputDir = "./out/"  
gridOutputFileName = "grid.Netcdf"  
geometryInputFileName = "geometry.toml"
```
## Memory allocation
Some parameters specifiying the allocation of the memory required for the grid generation have to be determined in theproperty file. Heap memory is specified by the property `scratchSize`. Based on experience for the grid generation `scratchSize = 5.0` is sufficient in most cases. Additionally, the [cell collector](#dvCollector) has to be allocated by defining the property `maxNoCells` which specifies the maximum number of cells per rank. Keep in mind that up to the `minLevel` no partitioning of the grid is performed. Therefore, every rank must be able to generate least all cells at the `minLevel`. We will elaborate further on the meaning of the various cell levels [below](#ugMeshRefinement). More details on the parallelization procedure are given [here](#nmParallelization). 

## Providing a geometry
Details on the files providing the actual geometry are given in the `geometry.toml`-file. The geometry of the 
computational domain is provided via ASCII .stl-files for a 3D setup and simple coordinate sequences in 2D. The paths to these 
files are given by setting the incremented property `filename.i` where i refers to the i-th file. For every new 
boundary condition applied to a part of the geometry a new file has to be provided. The total number of files is then 
given in `noSegments`. An example for a 3D domain the boundaries of which are given by three segments and which contains 
an additional cylinder which is given by a fourth file could look as follows:  
```
##### geometry.toml #####
filename.0 = "inflow.stl"  
filename.1 = "outflow.stl"  
filename.2 = "wall.stl"  
filename.3 = "cylinder.stl"  
noSegments = 4
```
The given segment files can be summarized to the the domain boundaries and, e.g, bodies within the computational domain by
listing the increments of the files in a respective field of the property `body_segments`:  
```
##### geometry.toml #####
body_segments.boundary = [0,1,2]  
body_segments.cylinder = [3]
```
@note Make sure your geometry is waterproof! If this is not the case, the grid generation algorithm cannot correctly determine whether a cell is inside or outside the fluid domain and the resulting grid will be corrupted.

Eventually, boundary conditions are assigned to each boundary defined by the segment files by setting `BC.i` to the ID 
of the desired boundary condition for the .stl with the respective increment i. For the example given above, this could 
look as follows:
```
##### geometry.toml #####
BC.0 = 1001    
BC.1 = 1002  
BC.2 = 100100  
BC.3 = 3003  
```
Now all geometries are defined and have been assigned a boundary condition. With the given information the 
`geometry.toml`-file is now complete.

## Mesh refinement {#ugMeshRefinement}
The mesh generation is based on the sequential subdivision of a base cube/rectangle (from now on simply refered to 
as cube) with the edge length \f$L_0\f$. The edge length \f$ L_0 \f$ is given by the largest dimension of the provide 
geometry segment files, such that the geometry fits entirely into the base cube. The cell length at a given level
\f$ l \f$ is then given by:  
\f{align}{
L_{l} = L_0 \left( \frac{1}{2} \right)^{l}
\f}

### Uniform refinement
The uniform refinement of the mesh is defined by the properties `minLevel` and `maxUniformRefinementLevel`.
-   `minLevel`: defines the refinement level up to which the grid generation runs serially. In the most simple case 
the partitioning of the grid takes place on the `minLevel` both in the solver and in the mesh generator. A high 
`minLevel` can be desirable to achieve a well balanced partitioning. The mesh generation of such grids can, however, 
require large amounts of memory as no parallelization is done before the `minLevel` is reached. More details on the parallelization are given [here](#nmParallelization).
-   `maxUniformRefinementLevel`: defines the level up to which the complete mesh is refined uniformely. Must be equal to or larger than the `minLevel`.

### Local refinement
Additional local refinement can be applied during mesh generation to properly resolve relevant parts of the 
computational domain. The maximum refinement level can be defined by setting the property `maxRfnmntLvl` accordingly. 
Adaptive mesh refinement strategies that are applied at runtime of the simulation are described [here](#nmAMR). 
For the purpose of a-priori mesh refinement two options are available:
-   boundary refinement
-   patch refinement

To enable the use of these refinement methods the property `localRfnMethod` offers the following options:
<table>
<tr><th> `localRfnMethod` </th>   <th> Description </th></tr>
<tr><td> <center> 0 </center> </td> <td> no local mesh refinement is applied </td></tr>
<tr><td> <center> 1 </center> </td> <td> only patch refinement is applied </td></tr>
<tr><td> <center> 2 </center> </td> <td> only boundary refinement is applied </td></tr>
<tr><td> <center> 3 </center> </td> <td> both patch and boundary refinement is applied </td></tr>
</table>

In the following we will go through the different options.

#### Patch refinement
Various options of differently shaped patches are available to refine the mesh within the patch region. The 
geometric shape of the patch is defined by providing the according token in the property `localRfnLvlMethods`. The 
required geometric paramters to define the shape of the patch are provided as comma-separated values in the property 
`localRfnLevelProperties`. An overview of the available patch geometries and the required parameters are given in in 
the following table:
<table>
<tr><th> Shape </th>    <th>`localRfnLvlMethods`</th>    <th>Parameters</th></tr>
<tr><td> Box </td>      <td><center> "B" </center> </td>      <td>\f$ x_{min},\, y_{min},\, z_{min},\, x_{max},\, y_{max},\, z_{max} \f$</td></tr>
<tr><td> Sphere </td>   <td><center> "R" </center> </td>      <td>\f$ x_{c},\, y_{c},\, z_{c},\, r\f$</td></tr>
<tr><td> Cylinder </td> <td><center> "C" </center> </td>      <td>\f$ x_{c,min},\, y_{c,min},\, z_{c,min},\ x_{c,max},\ y_{c,max},\, z_{c,max},\ r \f$</td></tr>
<tr><td> Tube </td>     <td><center> "T" </center> </td>      <td>\f$ x_{c,min},\, y_{c,min},\, z_{c,min},\ x_{c,max},\ y_{c,max},\, z_{c,max},\ r_{outer},\, r_{inner} \f$</td></tr>
<tr><td> Cone </td>     <td><center> "O" </center></td>      <td>\f$ x_{tip},\, y_{tip},\, z_{tip},\, x_{c,base},\, y_{c,base},\, z_{c,base},\, \varphi,\, r\f$</td></tr>
<tr><td> Hat (hollow cone) </td>  <td><center> "H" </center> </td>      <td>\f$ x_{tip},\, y_{tip},\, z_{tip},\, x_{c,base},\, y_{c,base},\, z_{c,base},\, \varphi,\, r,\, t \f$</td></tr>
<tr><td> Sliced cone (aligned)</td> <td><center> "A" </center> </td>      <td>\f$ x_{c,1},\, y_{c,1},\, z_{c,1},\, x_{c,2},\, y_{c,2},\, z_{c,2},\, r_1,\, r_2 \f$ </td></tr>
<tr><td> Sliced cone (non-aligned)</td> <td> <center> "N" </center></td>      <td>\f$ x_{c,1},\, y_{c,1},\, z_{c,1},\, x_{c,2},\, y_{c,2},\, z_{c,2},\, r_1,\, r_2,\, n_{1,x},\, n_{1,y},\, n_{1,z},\, n_{2,x},\, n_{2,y},\, n_{2,z} \f$ </td></tr>
<tr><td> Angled rectangular cuboid </td>  <td><center> "S" </center> </td>  <td>\f$ x_{c,min},\, y_{c,min},\, z_{c,min},\, x_{c,max},\, y_{c,max},\, z_{c,max},\, h,\, w,\, \f$</td></tr>
<tr><td> Cartesian cylinder segment </td> <td><center> "W" </center> </td>  <td>\f$ x_{c,min},\, y_{c,min},\, z_{c,min},\, r,\, l,\, dir,\, \varphi_{min},\, \varphi_{max} \f$</td></tr>
</table>
Refinement patches can be combined arbitrarily by simply listing their corresponding tokens in `localRfnLvlMethods`. Multiple patches are grouped according to their refinement level. The individual groups for the particular levels are separated by a hyphen. The corresponding geometric parameters are simply listed consecutively in `localRfnLevelProperties`. The following example will refine one box-shaped region and one sphere-shaped region by one refinement level. An additional box-region nested in the first box is refined once more. 
```
##### properties.toml #####
localRfnLvlMethods = ["BS-B"]  
localRfnLevelProperties = [0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 0.5, 0.1, 0.1, 0.1, 0.9, 0.9, 0.9]
```
@note 
- The refinement patches are applied consecutively starting at the lowest level. Every cell can only be refined by one level at a time. Thus, refinement levels can not be skipped. If regions of higher refinement levels are desired, they have to be refined to all intermediate levels before. 
- Level skips between cells and their neighbors should strictly be avoided. The size of the patches may need to be adjusted to ensure an appropriate transition between levels of neighboring cells.

#### Boundary refinement
The boundary refinement allows the refinement of the interfaces between the fluid domain and embedded bodies or domain boundaries. The interface is automatically refined up to the level specified by `maxBoundaryRfnLvl`. The control of the thickness of the cell layers at the individual levels is done through the property `localBndRfnMethod`. Three options are available:
<table>
<tr><th> `localBndRfnMethod` </th>   <th> Description </th></tr>
<tr><td> <center> 0 </center> </td>  <td> The thickness of the layer at the `maxBoundaryRfnLevel` is specified by the property `localMinBoundaryThreshold`. The thickness all of intermediate layers is given by the property `smoothDistance`. Both properties are specified in number of cells at the `maxUniformRefinementLevel`. Thus, individual layers can not be thinner than a single cell at the `maxUniformRefinementLevel` </td></tr>
<tr><td> <center> 1 </center> </td>  <td> The thickness of the layer at the `maxBoundaryRfnLevel` is specified by the property `localMinBoundaryThreshold` and is specified in fractions of a cell at the `maxUniformRefinementLevel`. The thickness of all intermediate layers is given by the property `smoothDistance` and is specified in number of cells at the respective level.</td></tr>
<tr><td> <center> 2 </center> </td> <td> The thickness of the layer at the `maxBoundaryRfnLevel` is specified by the property `localBndRfnDistance` and is given in non-dimensional length units. The thickness of all intermediate layers is given by the property `smoothDistance` and specified in number of cells at the respective level. </td></tr>
</table>

## Cut-off boundaries
During mesh generation so-called cut-off boundaries can be created. These are not defined by a geometry, instead the mesh is cropped according to specified parameters. The resulting boundaries are treated separately from the ghost cell approach that is applied to regular boundaries. Details on the numerical background of cut-off boundaries is given [here](#). To create a cut-off boundary the appropriate cut-off mode has to be chosen first be specifying the `cutOff`-property. The available options are given in the table below.  

<table>
<tr><th> `cutOff` </th>     <th> Description </th> </tr>
<tr><td> <center> 0 </center> </td>   <td> No cut-off is applied. </td></tr>
<tr><td> <center> 1 </center> </td>   <td> Cut-off is applied at the `minLevel`. </td></tr>
<tr><td> <center> 2 </center> </td>   <td> Cut-off is applied at all levels <= `minLevel`. </td></tr>
<tr><td> <center> 3 </center> </td>   <td> Cut-off is applied at all levels >= `minLevel`. </td></tr>
<tr><td> <center> 4 </center> </td>   <td> Cut-off is applied at all levels.</td></tr>
</table>

The geometric shape at which the cut-off is applied can be specified by the property `cutOffMethod`. The corresponding geometric limits are then given by the property `cutOffCoordinates`. Additionally, a number of cell layers that protrude these limits can be specified by the property 'cutOffNmbrLayers' depending on the chosen `cutOffMethod`. The available options for these properties are summarized in the following table.

<table>
<tr><th> Shape </th> <th> `cutOffMethod` </th>  <th> `cutOffCoordinates` </th> <th> `cutOffNmbrLayers` </th></tr>
<tr> <td> Plane cut-off </td> <td> <center> "P" </center> </td> <td> \f$ x,\, y,\, z,\, n_x,\, n_y,\, n_z \f$ </td> <td> requires \f$ 1 \f$ value.</td></tr>
<tr> <td> Box cut-off </td> <td> <center> "B" </center> </td> <td> \f$ x_{min},\, y_{min},\, z_{min},\, x_{max},\, y_{max},\, z_{max} \f$ </td> <td> requires \f$ 2\cdot nDim \f$ values.</td></tr>
<tr> <td> Inverse box cut-off </td> <td> <center> "iB" </center> </td> <td> \f$ x_{min},\, y_{min},\, z_{min},\, x_{max},\, y_{max},\, z_{max} \f$ </td> <td> requires \f$ 2\cdot nDim \f$ values.</td></tr>
<tr> <td> Cylinder segment on \f$x\f$-axis </td> <td> <center> "C" </center> </td> <td> \f$ y_c,\, z_c,\, \Delta\varphi,\, \varphi_c \f$ </td> <td> requires 1 value. </td></tr>
</table>
Multiple cut-off methods can be combined arbitrarily. The respective tokens in `cutOffMethod` are separated by a hyphen. The required numerical parameters in `cutOffCoordinates` and `cutOffNmbrLayers` are consecutevily listed separated by commas.

One might be tempted to not care much about the geometry files defining the outer domain boundaries when using cut-off boundaries. However, it is advisable not to use excessively large geometries even if the domain is trimmed down using a cut-off as the grid generator will fill the complete geometry with cells at the desired level *before* the cut-off is applied. Therefore, this practice can create a huge memory overhead that might exceed the limitations of the available hardware. 

## Mesh validity for Lattice-Boltzmann
Before the generated mesh is written to the output directory the grid generator will check the mesh on suitability for use with the [Lattice-Boltzmann](#nmLB) solver. This check will detect cells that have more than one level jump to their edge- or space-diagonal neighbors. If any such cells are detected the grid generator will issue a warning and the user is advised to check the generated grid carefully and if necessery adjust the employed refinement strategy. 
@note Although, the [Finite-Volume](#nmFV) solver generally allows the mesh to be used despite the warning, it indicates low mesh quality and the grid might need revision.   

## Advanced practices
In the following some advanced concepts are presented that allow to fine-tune the grid generation concept. These concepts can help to improve the performance of the grid generation, the subsequent simulation and enable the full potential of m-AIA's multi-solver approach. 
### Reduction factor
As we saw [before](#ugMeshRefinement) the cell size is restricted to subdivisions of the base cube with the edge length \f$ L_0 \f$ and therefore bound to the bounding box defined by the provided geometry. Sometimes, however, it is desirable to adjust the cell size more precisely. For that purpose the `reductionFactor`-property is introduced, which can take values between \f$ 1 \f$ and \f$ 2 \f$. The reduction factor is multiplied to the base length \f$ L_0 \f$ and thereby alters the edge length of the base cube. The cell size of all higher level cells is affected accordingly. 

### Multi-solver grids
To enable the full potential of m-AIA, namely the [coupling](#dvCoupler) of multiple [solvers](#dvSolver), it is possible to generate multi-solver grids by setting the property `multiSolverGrid = true`. The number of solvers is then given by the property `noSolvers`. Properties of the individual solvers \f$X\f$ are then addressed by adding the suffix `.X` to the corresponding property name. If properties are defined for specific solvers it is necessary to also prvide a default value with the suffix `.default` that will be used by all other solvers without a specified value.

### Initialize bodies from Level-Set function
Embedded bodies can also be initialized from a [level-set function](#mmLS) later during the simulation. Nevertheless, these bodies have to be specified in the geometry file. To tell m-AIA to ignore the specified bodies during the mesh generation set the property `GFieldInitFromSTL = true`. The exact bodies from the geometry file can then be specified by their boundary condition by listing the according boundary condition IDs in the property `bodyBndryCndIds`.

### Partition level shift
Fundamentally, the partitioning of the grid takes place at the `minLevel`. If the grid contains a wide range of refinement levels it can be difficult to achieve a well-balanced partitioning. The introduction of the so-called [partition level shift](#nmParallelization) allows to locally partition the grid on a higher refinement level. The threshold at which a partition level shift is triggered can be specified by setting the properties `partitonCellOffspringThreshold` and `partitionCellWorkloadThreshold`.

### Dynamic load balancing
Often, regions of mesh refinement are restricted to only small parts of the computational domain and refined cells, which are outside of the domain are deleted. In the course of the parallel grid generation this can lead to a high concentration of cells on only a few ranks. Enabling the [dynamic load balancing](#nmDLB) algorithm allows the repartitioning of the computational domain during the refinement process. The dynamic load balancing can be activated by setting `dynamicLoadBalancing = true`.
@bug Using dynamic load balancing in the generation of grids that have a partition level shift can result in corrupted grids.


