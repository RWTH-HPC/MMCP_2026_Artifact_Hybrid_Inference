MAIA is a multi-physics PDE solver framework with a focus on fluid dynamics
equations. It is developed at the Institute of Aerodynamics of the RWTH Aachen
University in Germany.

# Quickstart

1. Clone this repository 
   * Run `git clone git@git.rwth-aachen.de:aia/MAIA/Solver.git`
   * Or directly with its submodules `git clone git@git.rwth-aachen.de:aia/MAIA/Solver.git --recurse-submodules`
2. Run `./configure.py ? ?` in the top directory.
3. Run `make`
4. Execute a setup using `./maia #property_file`

# Build and view the documentation

1. Clone this repository
2. Run `/configure.py ? ?` in the top directory.
3. Run `make doc`. This will take a few minutes since the whole code is scanned using Doxygen.
    * Or run the fast compilation with `make doc-fast` to scan only the doc-folder for any changes.
4. Open the file doc/documentation.html in a web browser.

# Start developing

1. Create a private fork.
2. Implement/Improve/Fix something.
3. Create a merge request.
    1. Include a complete canary report with your merge request.
4. Wait for approval of your changes.


# Installing

## Prerequisites

Automatic configuration files are provided for:

* AIA
* JUWELS (Forschungszentrum Juelich)
* JURECA (Forschungszentrum Juelich)
* HAWK (HLRS Stuttgart)
* Claix (RWTH Aachen)

The following software packages are required for installation:

* Parallel NetCDF >= 1.9
* OpenMPI >= 4.0
* FFTW (with MPI support) >= 3.3.2

## Build

Building MAIA requires Python 3.x and CMake >= 2.8.

The following compilers are currently supported and tested:

* GCC >= 9.1
* Clang >= 9
* Intel >= 19.1
* PGI >= 19.2

Older versions might need some fixes.

## Contact

Questions/Problems regarding GIT/GITLAB ask **Ansgar** or **Julian**.
