# Tutorial 3 (Grid Generation) # {#ugTutorial3cartesian}
[TOC]

In this tutorial, you will learn the basic concept behind the parallel grid generator and in the process create your own grids with different types of refinement and cut-offs.

For help, queries and suggestions contact: t.wegmann@aia.rwth-aachen.de . 

## Theory 

See [Theory & Implementation Grid Generator (old)](http://ldap2.aia.rwth-aachen.de/mediawiki-1.22.1/index.php5/ZFS:Theory_%26_Implementation_Grid_Generator) and [Grid generation (old)](http://ldap2.aia.rwth-aachen.de/mediawiki-1.22.1/index.php5/ZFS:Grid_generation)!

## Tutorial
Preparation
* Compile maia in production, e.g. <code>./configure.py gnu production</code>
* Download files for this tutorial: [WS3.zip](https://git.rwth-aachen.de/aia/MAIA/Solver/-/wikis/uploads/39ff4aa3764286fbdc093158382c9ef2/WS3.zip).
* Optional: Take a single twelve node interactive session, see [this tutorial](@ref ugTutorial2cartesian).

This tutorial assumes that
* You downloaded the files in the link above.
* You are located in the folder WS3
* You use vim as your editor (you can use any editor you want naturally)
* The directory should contain the following files
~~~
  user@WS3/> ls
  geometry.toml InitialSetup.png out/ properties_grid.toml solution/ stl/
~~~

### 1) Geometry File 
The _geometry.toml_ file controls the source of the grid. It is used to choose how and which STL files the grid generation will use and set boundary conditions for the geometry.

* Open the file <code>geometry.toml</code>
~~~
  user@WS3/> vi geometry.toml
~~~
~~~
  noSegments = 6                // Number of stl files used.
 
  filename.default = ""
  BC.default       = 0
 
  filename.0 = "stl/inflow.stl"  // For each boundary segment, write which file it corresponds to.
  filename.1 = "stl/outflow.stl"
  filename.2 = "stl/top.stl"
  filename.3 = "stl/bottom.stl"
  filename.4 = "stl/left.stl"
  filename.5 = "stl/right.stl"
 
  "body_segments.boundary" = [0, 1, 2, 3, 4, 5] // Numbering of the boundary segments.
 
  BC.0 = 0 // For each segment, select a boundary condition id (in this case 0 as we only generate grids without using them to run the flow solver).
  BC.1 = 0 
  BC.2 = 0 
  BC.3 = 0
  BC.4 = 0
  BC.5 = 0
~~~

* As seen in the *geometry.toml* file, the STL files are located in the */stl* folder.
~~~
   user@WS3/> cd stl/
   user@WS3/stl/> ls
   bottom.stl  inflow.stl  left.stl  outflow.stl  right.stl  top.stl
~~~

Find more information about:
* [What are STL files ?](https://en.wikipedia.org/wiki/STL_%28file_format%29)
* [How to create your own simple 3D stls with ParaView! (old)](http://ldap2.aia.rwth-aachen.de/mediawiki-1.22.1/index.php5/Num:visualization#Create_your_stl_geometry)

### 2) Geometry Exercise 
* Create a unit sphere (center (0.0, 0.0, 0.0), radius 1.0) ASCII stl in paraview based on the information from the link above.
* Add the stl of the sphere to the geometry.toml file and assemble it as an object.
* Compare your answers with the solutions.

### 3) Grid Property File 
The *properties_grid.toml* file controls the grid generation process.
It is used to pinpoint global and local refinements of the grid as well as cut-offs.

* Open the file **properties_grid.toml**
~~~
  user@WS3/> vi properties_grid.toml
~~~
* Most important parts of interest are:
~~~
 // -------------- GRID AND REFINEMENT PROPERTIES --------------
 maxNoCells                = 75000    // maximum number of cell each process can have.
 minLevel                  = 4        // globally refine the grid to level 4 before parallel distribution occurs.
 maxUniformRefinementLevel = 5        // globally refine the grid to level 5 after parallel distribution occurs.
 maxRfnmntLvl              = 6        // Continue with special refinements (patches and/or boundary) until level 6 is reached (in this case 1 level of special refinement).
 localRfnMethod            = 3        // Set the special refinement type: 0 deactivated, 1 patch only, 2 boundary only, 3 both.
 
 // Patch only property
 localRfnLvlMethods = "B"                                  // Do a patch refinement of the Box-type on the first special refinement level.
 localRfnLevelProperties = [-1.0, -1.0, -1.0, 0.0, 0.0, 0.0] // Coordinate of the box patch (xmin, ymin, zmin, xmax, ymax, zmax). Beware of the length of the variable when you change a patch!
 
 // Boundary refinement only property
 smoothDistance = 1                                       // Cell distance per level the grid has to be smoothed.
 localMinBoundaryThreshold = 2                            // On the last level (in this case level 6), a distance of 2 cells of size maxUniformRefinementLevel (in this case 5) will be refined.
 localRfnBoundaryIds = 6                                  // Which boundary id to perform the boundary refinement on (Specified in the geometry.toml).
 maxBoundaryRfnLvl = 6                                    // Do boundary refinement until level 6 is reached.
 
 // Cut-Off properties
 cutOff = 0                                               // Activate or deactivate cut-offs and specify the cut-off level.
 cutOffMethod = "B"                                       // Specify cut-offs type (in this case a Box cut-Off).
 cutOffCoordinates = [-1.0, -1.0, -1.0, 0.0, 0.0, 0.0]      // Specify cut-offs coordinates and box-dimensions. Beware of the length of the variable when you change a patch type!
~~~

### 4) Grid Property Exercise 
* Change the file properties_grid.toml to:
   1. globally refine the grid until level 6 before parallel distribution.
   2. globally refine the grid until level 8 after parallel distribution.
   3. change the box patch at level 9 to a unit sphere patch.
   4. activate Cut-Off at all Levels
* Compare your answers with the solution files and the solution-image.

**Hints**: 
1. You may need to adapt the boundary refinement level as you increase other refinements!
2. You may need to increase the maximum number of cells as you refine the grid. (e.g to 500000)
3. More information about the Cut-Off options can be found [here (old)](http://ldap2.aia.rwth-aachen.de/mediawiki-1.22.1/index.php5/ZFS:Grid_generation#Using_cut-offs). Most importantly it should be noted that the Cut-Off patch specifies the remaining cells! This can lead to errors if no cells lie within the Cut-Off Box.
4. The error-message concerning the mesh validity for LBM may be disregarded in this tutorial.

### 5) You're on your own! 
The folder contains the files needed to create a grid around a cube stl using 4 processes and 75000 cells per process.
* The dimensions of the cube are (-2.0, -2.0, -2.0) to (2.0, 2.0, 2.0).
* The cube consist of 6 2D square STLs which are used as solid wall boundary conditions.
* The initial serial refinement goes up to level 4.
* The initial parallel refinement goes up to level 5.
* The maximum refinement is level 5.
* Cut-offs, patches and boundary refinements are deactivated.
* Running MAIA will give the following grid:
![tutorial_3_cube](figures/tutorial_3_cube.png)

* Solve the following exercises and visualize the results in paraview.
  1. Create a sphere stl in ParaView with coordinates (1.0, 1.0, 1.0) and radius 0.5.
  2. Modify a copy of the template to use the created sphere with boundary refinement      at level 6 and 7.
  3. Modify a copy of the template to create a box patch at (0.5, 0.5, 0.5) to (1.0, 1.0, 1.0) on level 6.
  4. Modify a copy of the template to create a box cut off from (0.0, 0.0, 0.0) to (0.5, 0.5, 0.5).
  5. Everything together, modify a copy of template with
    * globally refine the grid until level 5 before parallel distribution.
    * globally refine the grid until level 6 after parallel distribution.
    * a boundary refinement at level 7 and 8 of the sphere created earlier.
    * a box patch from (-2.0, -1.0, -1.0) to (0.0, 0.0, 0.0) on level 7.
    * a cylinder patch with coordinate (-1.0, -1.0, -1.0) to (-0.5, -0.5, -0.5) with radius 0.25 on level 8.
    * a box cut off from (-1.0, -1.0, -1.0) to (1.0, 1.0, 1.0).
* Compare your answers with the solutions.

**Hints:**  
  * Refer to the grid generation [user guide (old)](http://ldap2.aia.rwth-aachen.de/mediawiki-1.22.1/index.php5/ZFS:Grid_generation) in case you need some in-depth information.
  * In some cases, you may need to increases the maximum number of cells allowed per process and/or the number of processes.


