# Structured Tutorial 1 (Setup) # {#ugTutorial1structured}
[TOC]

This tutorial describes the steps to set up the structured solver in MAIA, which is necessary and recommended to run a basic simulation in the following tutorials. For more information on the configuration of the structured solver, please refer to the corresponding user guide section [configuration of the structured solver](@ref ugSTRCD).

@note
Unlike the Cartesian solver part of MAIA, the structured solver cannot generate the grid itself, i.e., no `stl` folder and no stl-files as for the Cartesian solver part of MAIA, please refer to the [Cartesian tutorials](@ref ugTutorial1cartesian)). Therefore, the grid must be created beforehand and then provided to the structured solver part of MAIA. For more details on the numerical method, please refer to the following sections on the [structured grid](@ref nmStructuredGrid) and the [strcutured finite volume solver](@ref nmSTRCD).

## Installation

Refer to the instructions given in the [installation](@ref ugInstallationGuide) section or [quickstart](@ref ugGettingStarted) and complete the following tasks:

1. Clone a recent version of MAIA from **GIT**.

2. Compile a version of MAIA with **production** settings using the **GNU** compiler.


## Folders and Files

In the following tutorial you will encounter the following folders and files:  

<table>
<tr><th> Folder / file name     </th>  <th> Applications </th></tr>
<tr><td> grid.hdf5   </td>  <td> Grid file in HDF5 format to be read in by the solver. The name needs to be specified in the property file (see next row of this table) </td></tr>
<tr><td> properties_???.toml   </td>  <td> Configuration file used to define properties for the **solver setup, restart, output, etc.** </td></tr>
<tr><td> /out                   </td>  <td> Output folder which contains the simulation results as HDF5 solution files `SomeTimestep.hdf5`, where `SomeTimestep` is the integer number of the time step at which the file was written out </td></tr>
<tr><td> /auxdata   </td>  <td> Folder for auxdata files, which contain spatially resolved wall properties, e.g., \f$c_f \f$ </td></tr>
<tr><td> /boxes   </td>  <td> Folder for box output files, i.e., HDF5 files which contain sub-volumes of the primitive (and possibly other) variables  </td></tr>
</table>


## Tutorial Setup

After having downloaded and installed (configured and compiled) MAIA, proceed with the following steps for setting up the tutorials.

1. **Download** the files of the tutorials for the strcutured solver part of MAIA from the given link and **unzip** with the command
   ~~~
   cd /some/folder/
   unzip ~/Downloads/tutorialNo1.zip
   ~~~

    <br>

2. In the top directory of each tutorial folder, **link** your MAIA executable to the tutorial directory, so you can save space and easily access it from within the tutorial. Repeat this for every tutorial:  
   ~~~
   cd tutorialNo1
   ln -s /home/someuser/Solver/src/maia ./maia
   ~~~

3. Test if MAIA is correctly linked by opening MAIA's help:  
   ~~~
   maia -h
   ~~~

@note
In the beginning you can leave the property files unchanged, as these tutorials are only intended to give an overview of the simulation workflow. Once you have gone through the simulation workflow once, you can play around with the property settings.

