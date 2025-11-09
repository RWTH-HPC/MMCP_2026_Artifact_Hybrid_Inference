First steps with MAIA: The DG block
==================================


Overview
--------

1.  Overview
2.  Introduction
3.  Getting MAIA and the testcases
4.  Building MAIA
5.  Setting up a simulation
6.  Running MAIA
7.  Viewing the results
8.  How to
9.  ZACC
10. Further information


Introduction
------------

MAIA is a multi-physics PDE solver framework with a focus on fluid dynamics
equations. It is developed at the Institute of Aerodynamics (AIA) of the RWTH
Aachen University.

This is a tutorial for the novice MAIA user with a focus on the DG block.
It is self-contained in the sense, that from downloading the source code up to
viewing the results of a simulation, all necessary steps are described.

Further information can be found on the [AIA wiki][10].
Additional references concerning specific aspects will be given throughout the
paragraphs and at the end of this document.

[10]: http://ldap2/mediawiki-1.22.1/index.php5/Main_Page


Getting MAIA and the testcases
-----------------------------

For the source code and the testcases, version control is performed with
[subversion][30].
They are stored in a so-called repository.
In order to access the MAIA subversion repository, you have to be logged in
at one of the AIA computers.

The MAIA subversion repository is located at the following url:

    http://svn/maia

Stable code is stored in a certain folder inside the repository, the
so-called trunk.
The latest trunk and testcases are located at:

    http://svn/maia/trunk
    http://svn/maia/testcases

Copies of trunk code or testcases which contain new or improved code are stored
in other folders inside the repository, so-called branches.
Once they are confirmed to be stable, they will be reintegrated into the trunk.
Branches of the trunk and corresponding testcases are located at:

    http://svn/maia/branches/BRANCHNAME
    http://svn/maia/branches/tc_BRANCHNAME

[30]: https://subversion.apache.org/

### Typical workflow

1.  Check out a working copy of the latest trunk.

        svn co http://svn/maia/trunk

    Checking out the trunk in your home directory is fine.
2.  Check out a working copy of the latest testcases.

        svn co http://svn/maia/testcases

    Make sure to check out testcases in a scratch folder as running the
    simulations later may produce a lot of data.
3.  Update to the latest version if there were changes by others since the
    checkout/last update.

        svn update path/to/working/copy
    
By the way: If you are inside a working copy or one of its subdirectories, the
url of the respective repository root (`http://svn/maia`) may
be abbreviated by `^`.
E.g., `http://svn/maia/trunk` becomes `^/trunk`.


Building MAIA
------------

### Typical workflow

1.  Go to working copy of MAIA.

        cd path/to/maia

2.  Configure the build system to use default compiler and build type
    (`gnu debug`).

        ./configure.py ? ?

3.  Compile and link MAIA [with `NP` processors].

        make [-j NP]

Steps 2 and 3 are worthwhile being discussed in a little more detail.

### Configuring the build system

You can also configure the `COMPILER` and the `BUILDTYPE`:

    ./configure.py COMPILER BUILDTYPE

Available compilers include:

*   `gnu`
*   `clang`
*   `pgi`
*   `intel`

Available build types include:

*   `debug`
*   `production`
*   `extreme`

Get a complete list of available compilers and build types:

    ./configure.py

Get a complete list of optional arguments:

    ./configure.py -h

Selected (optional) arguments:

*   Enable OpenMP extensions (if supported by compiler).

        ./configure.py -fopenmp COMPILER BUILDTYPE

*   Enable the MAIA debugging facilities.

        ./configure.py --enable-debug COMPILER BUILDTYPE

### Compiling and linking

Start compiling and linking:

    make

Start compiling with NP processors and linking:

    make -j NP

By the way: *Never ever* just enter `make -j`.
It will impair the work of others logged into the same frontend.
To compile with the optimum number of processors, add the following alias to
your `~/.bashrc` file:

    alias m='make -j $(python -c "import multiprocessing; 
    print(multiprocessing.cpu_count())")'

Changes in the `~/.bashrc` file take effect as you log in the next time.
Enter `m` into the command line instead of `make -j NP`.

### References

Further information is available here: 

*   [AIA wiki][40]
*   [trunk/README.md][41]

[40]: http://ldap2/mediawiki-1.22.1/index.php5/MAIA:User_Guide#Installation
[41]: http://svn/maia/trunk/README.md


Setting up a simulation
-----------------------

### The testcases directory

1.  Go to working copy of testcases.

        cd path/to/testcases

2.  Go to DG subdirectory.

        cd DG

3.  Go to folder of respective test.

        cd folder/of/simulation

### Setting the properties

MAIA uses NetCDF property files. A (raw) property file has the following outline:

    netcdf PROPERTYFILE {
      dimensions:
          [here are some constants to be used]
      variables:
          [here are some variable definitions used in this file]
      data:
          [here come the variable assignments]
    }

There are three types of property files in each testcase:

*   `geometry.cdl`
*   `properties_grid.cdl`
*   `properties_run.cdl` and/or `properties_restart.cdl`

### geometry.cdl

In this file, the geometry is defined and the boundary conditions are specified.

Selected properties:

*   `int BC.INDEX`
    
    This variable specifies the boundary condition id for segment no. `INDEX`.

*   `int body_segments.boundary()`
    
    This variable specifies the segment indices for which boundary conditions
    are to be defined.
    
*   `char filename.INDEX()`

    This variable specifies the location of the file where the geometry of
    segment no. `INDEX` is defined.
    `INDEX` ranges from `0` to `noSegments-1`.

*   `int noSegments`
    
    This variable defines the number of segments (straight lines in 2D, planes
    in 3D) that make up the surface of your geometry.

### properties_grid.cdl

In this file, the properties for grid generation are specified.

Selected properties:

*   `int flowSolver`

    This variable has to be set to `0` for grid generation.
    Otherwise, the flow solver will be called immediately after the grid
    generation, which is not supported.

*   `int gridGenerator`

    For grid generation, always set this variable to `1`.

*   `int gridGeneratorParallel`

    For grid generation, always set this variable to `1`.

*   `int gridGenPar`

    This uses the "new" parallel grid generator.
    For grid generation, always set this variable to `1`.

*   `int maxNoCells`

    This variable specifies the maximum number of cells that can be stored
    *on each MPI rank*. It has to be adjusted depending on the approximate
    expected grid size.
    If the value is chosen (much) too high, the generator will allocate too much
    memory and might run out of memory (and thus crash).
    If it is too small, MAIA will abort.

*   `int multiBlock`

    Always set to `1`.

*   `int noDomains`

    This variable is only relevant for the automatic execution of the tests. It
    does not affect the grid generation process.

Refinement properties:

*   `int initialGeometricalRfnLvl`

    Level to which the grid is initially refined to. At this level the domain
    decomposition takes place when running the solver in parallel later. That
    is, you need at least as many cells on this level as the number of MPI ranks
    you want to use. Depending on the case, you probably want at least 10-500
    of these cells per rank. In MAIA parlance, these cells at the initial level
    are also called *min cells*.

*   `int maxBoundaryRfnLvl`

    Refinement level of cells at boundaries.

*   `int maxGeometricalRfnLvl`

    Level to which the grid is refined homogeneously after having been refined
    to `initialGeometricalRfnLvl`.

*   `int maxRfnLvl`

    Maximum possible level of any cell contained in the grid.

`initialGeometricalRfnLvl <= maxGeomentricalRfnLvl <= maxBoundaryRfnLvl
<= maxRfnLvl`

Rule of thumb for refinement:

*   Set `maxGeometricalRfnLvl`, `maxBoundaryRfnLvl`, and `maxRfnLvl` to the same
    value, which depends on your spatial resolution requirements.
*   Set `initialGeometricalRfnLvl` to a sensible value depending on your
    parallelization needs, but at least one level lower than the maximum
    refinement level.

Further information is available here:

*   [AIA wiki][50]

[50]: http://ldap2/mediawiki-1.22.1/index.php5/MAIA:User_Guide#Problem_Setup

### properties_run.cdl

In this file, the properties for the actual simulation are specified.

Selected properties:

*   `int flowSolver`

    This variable has to be set to `1` (see above).

*   `int gridGenerator`

    This variable has to be set to `0` (see above).

*   `int multiBlock`

    Always set to `1`.

*   `int noDomains`

    This variable is only relevant for the automatic execution of the tests. It
    does not affect the solver run.

Solution properties:
    
*   `int analysisInterval`

    This variable denotes the time step interval after which a summary is
    printed on the screen.
    
*   `int restartInterval`

    This variable denotes the time step interval after which a restart file is
    written (see 'Continue simulation from restart file' below).

*   `int solutionInterval`

    This variable denotes the time step interval after which a solution file is
    written.

*   `int timeSteps`

    This variable denotes the maximum number of timesteps. When reached, the
    simulation exits gracefully.

Spatial discretization:
    
*   `char blocktype()`

    For DG simulations, this variable is set to `MAIA_DISCONTINUOUS_GALERKIN`.

*   `char dgIntegrationMethod()`

    This variable is set to `DG_INTEGRATE_GAUSS` or
    `DG_INTEGRATE_GAUSS_LOBATTO` and denotes the used integration scheme. If
    there is no specific reason to do it differently, always use
    `DG_INTEGRATE_GAUSS`.

*   `int initPolyDeg`

    This variable specifies the used polynomial degree.

*   `int spaceDimensions`

    This variable is set to `2` or `3` for 2D or 3D simulations respectively.
    It should be set to the same value as `nDim` in the `dimensions` section of
    the property file.

Time integration:

*   `double cfl`

    This variable specifies the CFL number.
    As a rule of thumb, set it to `1.0`.
    If your simulation crashes or produces `NaN` ('not a number') values, it
    might be because your timestep size exceeded the stability limit.
    In this case, try lowering the CFL number.

*   `double finalTime`

    This variable specifies the finishing time.
    `finalTime` along with `cfl` determine the number of timesteps of your
    simulation.
    When reached, the simulation exits gracefully.

*   `double startTime`

    This variable specifies the starting time.
    Typically set to `0.0`.

Flow variables and initialization:

*   `char dgSystemEquations()`

    This variable specifies the set of quations.
    For acoustic perturbation equations, set to `DG_SYSEQN_ACOUSTICPERTURB`.
    For linear scalar advection equations, set to `DG_SYSEQN_LINEARSCALARADV`.

*   `int initialCondition`

    This variable specifies the initial condition for the respective set of
    equations.
    Set to one of the cases for [acoustic perturbation equations][51] or [linear
    scalar advection equations][52].

[51]: http://svn/maia/trunk/src/maiadgsyseqnacousticperturb.h
[52]: http://svn/maia/trunk/src/maiadgsyseqnlinearscalaradv.h

Grid and refinement properties:

*   `int maxNoCells`

    Maximum number of Cartesian cells (see above).
    This value is meant per processor.

*   `int maxNoSurfaces`

    Maximum number of surfaces between cells.
    This value is meant per processor.
    As a rule of thumb, set it to `nDim` times `maxNoCells`.

Further information is available here:

*   [AIA wiki][53]

[53]: http://ldap2/mediawiki-1.22.1/index.php5/MAIA:User_Guide#Solver_run


Running MAIA
-----------

### Using the AIA cluster

At the AIA, simulations should not be run on the personal computers but on the
AIA cluster. The cluster consists of the frontends (login nodes named `fe1` to
`fe6`), from which the actual compute nodes can be reached. There are 3 kinds
of compute nodes:

*   Nodes with 8 cores (queue `eight`)
*   Nodes with 12 cores (queue `twelve`)
*   Nodes with 24 cores (queue `twentyfour`)

On each node, there are 2 GB of memory available per core.
If a simulation requires more than that, the cores have to be distributed
across more nodes.
The processors on the 12-core nodes allow hyperthreading, i.e., you can have up
to 24 virtual cores running on the 12 physical cores. The 24-core nodes also
have hyperthreading.

There are three queues through which each kind of node can be accessed:
`eight`, `twelve` and `twentyfour` .
At the AIA, we use [TORQUE][60] as our resource manager.
Instead of using `qsub`, there are also job scripts available for interactive
sessions.
Execute

    qi NOCORES [NONODES]

where the number of cores denoted by `NOCORES` is either `8`, `12` or `24`,
depending on the desired queue.
Thus, the command is either `qi 8`, `qi 12` or `qi 24`.
If `NONODES` is not specified, its default value is `1`.
To activate X forwarding, add the `-X` flag (e.g. `qi -X 8`).

[60]: http://www.adaptivecomputing.com/products/open-source/torque/

Thus, an interactive session may be started like this:

1.  Check for available nodes on the cluster.

        pbshosts

2.  Start an interactive session with 3 nodes at 12 cores each.

        qi 12 3

Further information is available here:

*   [AIA wiki][61]

[61]: http://ldap2/mediawiki-1.22.1/index.php5/Com:cluster

Once an interactive session is started, you can either run MAIA manually or use
canary.

### Running MAIA manually

Example workflow:

1.  Go to simulation folder.

        cd path/to/testcases
        cd DG
        cd folder/of/simulation

2.  Generate properties file for the geometry.

        ncmpigen -o geometry geometry.cdl

3.  Generate properties file for the grid.

        ncmpigen -o properties properties_grid.cdl

4.  Generate grid.

        mpirun -np NP path/to/maia/src/maia

    where `path/to/maia/` is the location of your MAIA working copy.
    `NP` is the number of MPI ranks you want to use.

5.  Generate properties file for simulation.

        ncmpigen -o properties properties_run.cdl

6.  Start simulation.

        mpirun -np NP path/to/maia/src/maia

    where `path/to/maia` is the location of your MAIA working copy.
    `NP` is the number of MPI ranks you want to use.

By the way, you can export an enviroment variable to save you the work of
specifying the executable every time you start a simulation:

    export FOO=path/to/maia/src/maia

Now, you start a simulation with:

    mpirun -np NP $FOO

Further information is available here:

*   [AIA wiki][62]

[62]: http://ldap2/mediawiki-1.22.1/index.php5/MAIA:Readme#Getting_started

### Using canary

canary is a testing infrastructure of MAIA at the AIA.
It lets you easily run the testcases at the AIA cluster via command line.

Example workflow:

1.  Go to working copy of testcases or one of its subdirectories.

        cd path/to/testcases/[subdirectory/]

2.  Start canary with MAIA executable.

        cry test -z path/to/maia/src/maia
	
Note that canary is a recursive script which means that all testcases in the
current directory or below will be executed.
If no executable is specified, canary will fetch it from the path stored in the
environment variable `MAIA`.

Further information is available here:

*   [AIA wiki][63]
*   [cry README][64]

[63]: http://ldap2/mediawiki-1.22.1/index.php5/MAIA:Canary:UserGuide
[64]: /home/jhe/canary/cur/README.md

### Using qcom
If no nodes on the cluster are free you can also use qcom as a simple and quick
alternative to using a interactive session and writing a full blown job script.
The command invocation is similar to `qi`:

    qcom NOCORES [NONODES] CMD

Just as with `qi` you choose the nodes with NOCORES and NONODES. CMD is the
command you want to execute on the nodes (e.g. `mpirun -n 24 maia`).
You can also pipe your commands into qcom:

    cat commands | qcom 12 3

But at this point you should probably write a proper job script.
Further information is available here:

*   [AIA wiki][65]

[65]: http://ldap2/mediawiki-1.22.1/index.php5/Com:cluster

Viewing the results
-------------------

### Paraview

The solution of your simulations can be visualized with [Paraview][70].
Paraview is an open-source, multi-platform data-analysis and visualization
application.
At the AIA, you can use it to visualize NetCDF files, i.e. grid, solution and
restart files.

[70]: http://www.paraview.org/

How to open a solution file:

1. Start Paraview with 

        paraview --mpi

2.  Open your solution file by selecting it and clicking 'OK'.

    'File' (in the menu bar) -> 'Open...' -> ...

3.  On the left side, there should be a tab called 'Properties'.
    Press the 'Apply' button.

4.  In the 'Properties' tab, there should be a section called 'Coloring'.
    Pick the variable you want to have displayed.

How to plot over a line:

1.  Select the solution file in the pipeline browser and create a
    'Plot Over line' filter.

    'Filters' (in the menu bar) -> 'Data analysis' -> 'Plot Over Line'

2.  In the 'Properties' tab, enter the starting point ('Point1') and endpoint
    ('Point2') coordinates of your desired line.
    Then, press the 'Apply' button.

3.  In the 'Properties' tab, there should be a section called 'Series
    Parameters'.
    Select or unselect according to the variables you want to have displayed.

How to extract a plane:

1.  Select the solution file in the pipeline browser and create a
    'Slice' filter.

    'Filters' (in the menu bar) -> 'Common' -> 'Slice'

2.  In the 'Properties' tab, enter the coodinates of the origin and the normal
    vector.
    Then, press the 'Apply' button.

3.  In the 'Properties' tab, there should be a section called 'Coloring'.
    Pick the variable you want to have displayed.

### Gnuplot

In MAIA, you have the option to generate additional output in the form of line
data (see "How to" below).
If there is line data, you can visualize it with [gnuplot][71].
For example, plot the sixth column over the first:

[71]: http://www.gnuplot.info/

1.  Start gnuplot.

        gnuplot
    
    Make sure that X forwarding is activated for your SSH session (see above),
    i.e., activate the `-X` option when using the `ssh` command or the `qi`
    command.

2.  Plot the sixth column over the first for a file located at
    `path/to/file`.

        plot "path/to/file" using 1:6 with lines

3.  Leave gnuplot.

        quit

For a 2D simulation of the acoustic perturbation equations:

*   Columns 1 to 2 represent the x, y coordinates respectively.
*   Columns 3 to 4 represent the u, v velocity components respectively.
*   Column 5 represents the pressure.
*   Column 6 represents `spongeEta` if `writeSpongeEta` is set to `1`.

For a 3D simulation of the acoustic perturbation equations:

*   Columns 1 to 3 represent the x, y, z coordinates respectively.
*   Columns 4 to 6 represent the u, v, w velocity components respectively.
*   Column 7 represents the pressure.
*   Column 8 represents `spongeEta` if `writeSpongeEta` is set to `1`.

For a 2D simulation of the linear scalar advection equations:

*   Columns 1 to 2 represent the x, y coordinates respectively.
*   Columns 3 represents the scalar variable.

For a 3D simulation of the linear scalar advection equations:

*   Columns 1 to 3 represent the x, y, z coordinates respectively.
*   Columns 4 represents the scalar variable.

By the way:
While gnuplot offers interactive sessions, it is also possible to store commands
in a script file and execute them in batch mode.
This is highly recommended for anything more complex than a one-liner.
In order to execute the commands stored in the script file, enter the following
line (into the terminal):

    gnuplot path/to/script/file


ZACC
----

ZACC (MAIA Automatic Checking Collection) consists of several scripts that allow
you to automatically compile and run MAIA with the testcases.
It lets you check several configurations at once, i.e. different compilers,
different optimization levels and Valgrind enabled/disabled.

For each trunk commit, the corresponding revision will be tested against
several configurations.
Afterwards, a report is sent to all MAIA users with the current testcase status
and if there was a change in the results since the previous automatic check.

Thus before you check in new code, you want to verify it against ZACC manually:

    zacc_run.py

The following command line arguments are worthwhile being specified:

*   --svn-url-maia
    This is the url of your MAIA branch.
*   --svn-url-testcases
    This is the url of your testcases branch.
*   --mail-to
    This is the email address where the ZACC report is sent to.

Good advice:

*   Make sure that you execute `zacc_run.py` from within a scratch folder.
*   Check how busy the AIA cluster is (`pbshosts`) and how long the queue is
    (`qstat`).
*   Run ZACC only during the night or on weekends.
    (Unless you only run certain setups, see `zacc_run.py -h`.)
*   Consider to use [screen][80] (see "How to" below).
*   Have a look at the [AIA wiki][81].

[80]: http://www.gnu.org/software/screen/
[81]: http://ldap2/mediawiki-1.22.1/index.php5/MAIA:Testcases#MAIA_Automatic_Checking_Collection_.28ZACC.29


How to
------

### Continue simulation from restart file

Similar to solution files, a restart file stores the information on the flow
field at a certain timestep.
It can be used to start a simulation not from the initialized state, but at an
intermediately converged state.
This might come in handy, e.g. when the simulation has to be interrupted or
when boundary conditions have to be exchanged.

Example workflow:

1.  Set `restartFile` in `properties_run.cdl` to `1`.

2.  Set `restartTimeStep` in `properties_run.cdl` to the correct number of
    timesteps already calculated.
    Note that there has to be a restart file named `restart_TIMESTEPS.Netcdf`
    in the `out/` folder where `TIMESTEPS` denotes the number of timesteps
    alread calculated.

3.  Run simulation as usual.

### Write out line data

If you are not interested in the entire flow field, line data might be a viable
option for you.
A line is specified by a starting point, an endpoint and a number of intervals.
If this feature is activated, the state of the flow field at each interval over
the line is written to a file.
Thus, you can use gnuplot to visualize the data from the file in a classic x-y
diagram (see above).

Declare and set the following variables in `properties_run.cdl`:

*   `int saveLineDataFile`

    Set to `1`.

For each line indexed by `INDEX` (count starts from 0):

*   `double endPoint_INDEX()`

    This variable contains the comma-separated coordinates of the endpoint of
    line `INDEX`.

*   `char lineDataFileName_INDEX()`

    This variable specifies the name of the file where the line data is to
    be written.

*   `int noIntervals_INDEX`

    This variable denotes the number of intervals on line `INDEX`.

*   `double startPoint_INDEX()`

    This variable contains the comma separated coordinates of the starting
    point of line `INDEX`.

Note that your lines must not coincide with one of the coordinate axes.

An example of line data is given by the
[2D acoustic perturbation equation testcase for a pressurepulse][90].

[90]: http://svn/maia/testcases/DG/2D_square_acousticperturb_pressurepulse_np4_linedata/properties_run.cdl

### Use sponge layers

In external flow simulations, exact boundary conditions are usually unknown.
In order to artificially suppress unphysical reflections at the external
boundaries, a sponge layer is used.
It consists of an area close to the external boundaries where a negative source
term on the right-hand side becomes active that forces the state variables
towards a target value.

In MAIA, you can assign a sponge layer with its respective thickness to each
boundary condition individually.
`spongeEta` is `0` in the non-sponge area and `1` at the external boundary.
In between, `spongeEta` is interpolated by a simple square function.

The sponge feature is *supported only by the acoustic perturbation equations*.

Declare and set the following variables in `properties_run.cdl`:

*   `double defaultSpongeThickness`

    Sponge thickness for all boundary conditions listed in `spongeBCs` without
    `spongeThickness_<id>`.

*   `int spongeBCs()`

    Comma separated list of boundary condition ids.
    Only for boundary conditions listed in `spongeBCs`, a sponge layer is
    created.
    For each sponge layer, a thickness has to be specified (see
    `spongeThickness_ID`).
    If none is specified, a default thickness has to be specified (see
    `defaultSpongeThickness`).

*   `double spongePressureInfy`

    Set the desired farfield pressure for the acoustic perturbation equation.

*   `double spongeSigma`

    `spongeSigma` represents a damping coefficient.
    It has the same value for all boundary conditions.
    A sensible value would be between `0.3` and `1.0`.
    If none is specified, a hard-coded default value of `1.0` will be used.

*   `double spongeThickness_ID`

    Sponge thickness for boundary condition identified by `ID`. 

*   `int useSponge`

    Set to `1` to activate sponge mechanism.

If you want `spongeEta` to be written into the solution files:

*   `int writeSpongeEta`

    Set to `1`.

A sponge example is given by the
[2D acoustic perturbation equation sponge testcase][91].

[91]: http://svn/maia/testcases/DG/2D_square_acousticperturb_sponge_np8/properties_run.cdl

### Format source code

Selected formatting with vim:

1.  Edit your `~/.vimrc` file and enter the following lines:

        " Commands to invoke clang-format "
        map <C-K> :pyf /home/mic/.pool/.src/llvm/tools/clang/tools/clang-format/clang-format.py<CR>
        imap <C-K> <ESC>:pyf /home/mic/.pool/.src/llvm/tools/clang/tools/clang-format/clang-format.py<CR>

2.  Open source code file in vim.

3.  Mark the parts of the code to be formatted and press `CTRL-K`.

Formatting a whole file:

1.  Use clang-format and specify location of source code file.

        clang-format -i path/to/file --style=file

Do not format legacy code as this will only cause noise the next time you 
submit your changes.
Only format lines that you changed anyways.

General advice:

*   Read style guide in [AIA wiki][92].
*   Use other code as reference.

[92]: http://ldap2/mediawiki-1.22.1/index.php5/MAIA:Development:Guidelines

### Start long-running simulations

If want to start long-running simulations during an interactive session, you
normally have to stay logged in for the whole duration.
Use [screen][93] to leave and later resume your session without interrupting
your simulations in the meantime.

[93]: https://www.gnu.org/software/screen/

Screen is a full-screen window manager that multiplexes a physical terminal
between several processes, typically interactive shells.
Here, a typical workflow is given:

1.  Invoke screen [with a session name of your choice].

        screen [-S sessionname]

2.  Start your simulations.

3.  Leave screen session to be resumed later.
    Programs continue to run in the meantime.
    Press the following keys during your screen session:

        CTRL+a CTRL+d

4.  Print a list of session identification strings.

        screen -ls

5.  Pick the identification string of your screen session from the printed list
    and resume your screen session.

        screen -r session.identification.string

6.  Repeat steps 2 to 5 as often as needed.

7.  Kill the screen session.
    Press the following keys during your screen session:
        
	      CTRL+a CTRL+k

    Alternatively, logging out (using `CTLR+d` or `exit`) will also terminate
    the screen session.

Further information is available here: 

*   [GNU screen manual][94]

[94]: http://www.gnu.org/software/screen/manual/


Further information
-------------------

Items marked by (**) are available from Michael or Vitali.

### General

*   [trunk readme][1000]
*   [testcases readme][1001]

[1000]: http://svn/maia/trunk/README.md
[1001]: http://svn/maia/testcases/README.md

### On using Linux

*   [bash][1010]
*   [vim][1011]
*   [svn][1012]

[1010]: http://www.gnu.org/software/bash/manual/
[1011]: http://vimdoc.sourceforge.net/
[1012]: http://svnbook.red-bean.com/

### On C++

*   [A wiki book][1020]
*   [C and C++ reference][1021]
*   [Google C++ style guide][1022]

[1020]: http://en.wikibooks.org/wiki/C%2B%2B_Programming
[1021]: http://en.cppreference.com/w/
[1022]: https://google-styleguide.googlecode.com/svn/trunk/cppguide.html

### On parallelization via MPI/OpenMP

*   [Presentations by FZ Jülich][1030]
*   [Tutorial material by Argonne National Laboratory][1031]
*   [Training materials by Lawrence Livermore National Laboratory][1032]
*   [MPICH guides][1033]
*   [MPICH manpages][1034]
*   [OpenMPI documents][1035]
*   [MPI documents][1036] (including standards)

[1030]: http://www.fz-juelich.de/ias/jsc/EN/Expertise/Services/Documentation/presentations/_node.html
[1031]: http://www.mcs.anl.gov/research/projects/mpi/tutorial/
[1032]: https://computing.llnl.gov/?set=training&page=index
[1033]: http://www.mpich.org/documentation/guides/
[1034]: http://www.mpich.org/documentation/manpages/
[1035]: http://www.open-mpi.org/doc/
[1036]: http://www.mpi-forum.org/docs/docs.html

### On the DG method

*   Gassner: Discontinuous-Galerkin-Verfahren. 2013 (**)
*   Schlottke/Cheng/Lintermann/Meinke/Schröder: A direct-hybrid method for computational aeroacoustics. AIAA Conference, Dallas, June 22-26, 2015 (submitted) (**)
*   Hindenlang/Gassner/Altmann/Beck/Staudenmeier/Munz: Explicit discontinuous Galerkin methods for unsteady problems. Computers & Fluids, 2012, 61, 86-93 (**)
*   Kopriva: Implementing Spectral Methods for Partial Differential Equations. Springer, 2009 (**)
*   Hesthaven/Warburton: Nodal Discontinuous Galerkin Methods. Springer, 2008 (**)

### On visualization with ParaView/Gnuplot

*   [Official gnuplot documentation][1050]
*   [Gnuplot In Action][1051] (**)
*   [Paraview: User guide, documentation, tutorials, etc.][1052]

[1050]: http://www.gnuplot.info/documentation.html
[1051]: http://www.manning.com/janert/
[1052]: http://www.paraview.org/resources/
