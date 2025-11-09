MAIA with LIKWID
===============

Introduction
------------
LIKWID (Like I Know What I am Doing) is a performance tool that works on Intel
and AMD processors on the Linux operating system[[1]]. It can either be used as
a wrapper to measure the execution of an application itself or with its Marker
API, which allows the profiling of single parts of the code. The raw hardware
measurements are interpreted by so-called *performance groups*, which are
configurations that contain the instructions for different processors.

Implementation
--------------
The performance measurement of MAIA with LIKWID is implemented by using its
Marker API together with the MAIA timer class. This allows the evaluation of
Flops for specific parts of the code. The start and stop routines of the LIKWID
counters are added to the existing timer functions, with the unique
identification string being obtained from the already assigned timer id. The
actual measurement is performed by the tool `likwid-perfctr`, which evaluates
all counters and stores their results in a file. These values are afterwards
inserted into the MAIA log to obtain a single output file. At the current state
of the development the profiling can only be performed for single-threaded
executions.

Prerequisites
-------------
The main requirement is the installation of LIKWID, including its access daemon.
A detailed build instruction can be found on the LIKWID
[wiki](https://github.com/RRZE-HPC/likwid/wiki). Additionally the path to
LIKWID needs to be added to the `PATH` environment variable:

    PATH=$PATH:/pds/opt/likwid/bin
    export PATH

For the AIA, all necessary information on likwid was added to the CMake build
system and the likwid executables are already on the `PATH`.

Examples
--------
The measuring of double precision Flops and the merging of their log files is
handled by the *measureFlops.py* script. It is located in `auxiliary/` and
expects the path to the MAIA executable as an argument. For example, to start a
MAIA simulation with a concurrent evaluation of the Flops, execute:

    measureFlops.py ../maia_likwid/maia

This yields the log file `maia_log_flops`, which includes the measured flops next
to the corresponding timer results.

    [100.0%] Total                   2.439828 [sec]   1.147660e-05 [DP MFlops/s]
       [98.5%] Run Environment       2.403192 [sec]   2.496737e-05 [DP MFlops/s]
       ...

Please note that parallel execution of MAIA with likwid is not supported.

[1]: https://github.com/RRZE-HPC/likwid
