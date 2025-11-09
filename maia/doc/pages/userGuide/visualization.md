# Visualization # {#ugVisualization}

For the visualization of restart and solution files the open source software ParaView by Kitware is used, see https://www.paraview.org. To visualize your data, you, therfore, need to install ParaView.
When installing ParaView, make sure to enable MPI-, Python- and QT-support.

## Solution files (VTU)
The @ref ugFV allows to write solution files in VTU format by defining the property `outputFormat = "VTU"` in the properties file. VTU is a ParaView native file format, therefore, VTU files can be loaded by any standart ParaView installation. Just click load, select the desired solution file and hit apply.

## Restart files (Netcdf)
In @maia the restart files and the grid file are written in Netcdf format. To load Netcdf files created with @maia, the @maia ParaView plugins from the git repository are required.
Several custom @maia plugins are available.
The most important is the @maia reader plugin called AIAReaderParrotEmbedded. The reader plugin is an interface between the @maia file format and ParaView. It reads the @maia file and converts the Cartesian mesh to a vtkunstructuredgrid, which is then passed to ParaView.
The other plugins are independet of @maia and present custom ParaView filters usefull for flow visualization. A list of the available plugins is given below.

## Parallel visualization
The @maia reader plugin is fully parallelized. To load data files in parallel, the ParaView Client musst be connected to a ParaView Server.
For the parallel visualization the @maia reader plugin adds halos cells to the grid for the communication between the grid boundaries. The additional memory required for the halo cells is computed automatically.
However, if a large number of processors is used, a manual adjustment of the precalculated memory increase is necessary. Therefore, if your ParaView session crashes, you can try increasing the property
```
Increase average (cell) memory by: X.X
```
in the ParaView GUI.

##### List of available plugins.
<table>
<tr><th>Plugin                      </th><th>Type  </th> <th>Description </th></tr>
<tr><td>AIAReaderParrotEmbedded     </td><td>Reader</td> <td>The reader plugin is used to load Netcdf files created with  </td></tr>
<tr><td>Calc                        </td><td>Filter</td> <td>Transfers cell data to point data and stores velocity components u, v, w in the vector velocity </td></tr>
<tr><td>@f$\Delta@f$-Criterion      </td><td>Filter</td> <td>Computes delta criterion
\f$ \Delta = (\frac{1}{3} Q)^3 + (\frac{det \nabla u}{2})^2\f$,
see Q-Criterion for \f$ Q\f$
</td></tr>
<tr><td>Partition                   </td><td>Filter</td> <td>Visualizes the grid partitioning for a given number of processes. Can only be used on grid files </td></tr>
<tr><td>Kinematic Vorticity         </td><td>Filter</td> <td>Computes the kinematic vorticity </td></tr>
<tr><td>@f$\lambda_2@f$-Criterion   </td><td>Filter</td> <td>Computes the @f$\lambda_2@f$-criterion  
(Second eigenvalue \f$ \lambda_2 \f$ of (\f$\Omega^2 + S^2 \f$), see Q-Criterion for \f$ \Omega \f$ and \f$ S \f$)
</td></tr>
<tr><td>Q-Criterion                 </td><td>Filter</td> <td>Computes the Q-criterion
\f$ Q = \frac{1}{2} \cdot (||\Omega ||_2^2 - ||S ||_2^2) \f$  
\f$ \Omega = [\partial u_i / \partial x_j - \partial u_j / \partial x_i] \, ,
\quad S = [\partial u_i / \partial x_j + \partial u_j / \partial x_i] \f$
</td></tr>
<tr><td>Solution Diff               </td><td>Filter</td> <td>Compares to data files </td></tr>
<tr><td>Heat Release                </td><td>Filter</td> <td>Deprecated </td></tr>
<tr><td>Surface Error Calculator    </td><td>Filter</td> <td>Deprecated </td></tr>
<tr><td>Surface To Grid             </td><td>Filter</td> <td>Deprecated </td></tr>
<tr><td>Triangulate Boundary Filter </td><td>Filter</td> <td>Deprecated </td></tr> 
<tr><td>WSS                         </td><td>Filter</td> <td>Deprecated </td></tr>
</table>


##### Building the plugins

TODO: Move to plugin documentation

To build the reader plugin you need a working version of ParaView and the @maia source files. The necessary @maia io routines required to read and load the @maia Netcdf files are then directly compiled into the reader plugin.
First, check out the latest plugins from the git repository. To build the plugins use the **configure** script located in the main directory.
1. Create directories **build** and and **install**
2. Got to directory **build**
3. Update paths in the respective host file in <b>./aux/hosts/</b>  
   If you have not yet created a host file, you need add one for your machine and update the **aux/gethost** script
4. Run <b>../configure</b> from build directory
5. If an error occurs, check the config.log and solve the issue
6. Run <b>ccmake ../.</b> and update the *CMAKE_INSTALL_PREFIX* to directory **install**
7. Run **make install**

By running **make install** builds dynamic libraries <em>.so</em> and stores them in **install/lib64/**. You can now either load these libraries in ParaView manually or set the *$PV_PLUGIN_PATH* such that ParaView finds the libraries. Once the libraries are loaded, you can use them in ParaView.