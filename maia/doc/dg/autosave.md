End autosave feature in the DG block
================================


Introduction
-----------

The end autosave feature saves a restart file just before your reserved wall-time
runs out. This way you can restart your calculations right after they were
interrupted. Note that this feature currently only works in the DG block.


Usage
-----

To use the end autosave feature the environmental variable `MAIA_JOB_END_TIME` must
declare when the job ends. The format is of the variable is in unix
time(seconds since 1970-1-1 00:00:00 UTC).

To ease the user with the calculation of the end time a helper script is
supplied: You can find it under aux/enable_autosave.sh. To use it source it in
your job script as early as possible.  It will look for a wall-time duration in
the hh:mm:ss format in the script your sourcing it from. If no duration can be
found the user must specify it manually either in the hh:mm, hh:mm:ss or ssss
format.

There is one configuration property which can be set in the properties file:
`noMinutesEndAutoSave` specifies how many minutes before the job ends a restart
file will be written. If it is not specified the default value of 10 minutes is
used.

Note that on JUQUEEN the environmental variables are not automatically passed
to MAIA but must be specified manually:

    runjob --exp-env MAIA_JOB_END_TIME --exe path/to/maia

Otherwise this feature won't work.


Example job script for Hazel Hen
--------------------------------

    #!/bin/bash
    #PBS -N jobName
    #PBS -l nodes=1:ppn=24
    #PBS -l walltime=00:10:00

    cd $PBS_O_WORKDIR

    . maia/aux/enable_autosave.sh

    /opt/cray/netcdf/4.3.2/bin/ncgen -o geometry geometry.cdl
    /opt/cray/netcdf/4.3.2/bin/ncgen -o properties properties_grid.cdl
    aprun -n 24 maia/src/maia

    /opt/cray/netcdf/4.3.2/bin/ncgen -o geometry geometry.cdl
    /opt/cray/netcdf/4.3.2/bin/ncgen -o properties properties_run.cdl

    aprun -j 2 -n 24 -N 24 maia/src/maia

Example job script for JUQUEEN
------------------------------

    #@ shell = /bin/bash
    #@ job_name = jobName
    #@ output = $(job_name).$(jobid)
    #@ error = $(output)
    #@ environment = COPY_ALL
    #@ job_type = bluegene
    #@ bg_size = 32
    #@ wall_clock_limit = 0:10:00
    #@ queue

    . maia/aux/enable_autosave.sh

    /bgsys/local/netcdf/bin/ncgen -o geometry geometry.cdl
    /bgsys/local/netcdf/bin/ncgen -o properties properties_grid.cdl
    runjob --np 32 --exe maia/src/maia

    /bgsys/local/netcdf/bin/ncgen -o geometry geometry.cdl
    /bgsys/local/netcdf/bin/ncgen -o properties properties_run.cdl

    runjob  --ranks-per-node 16 --exp-env MAIA_JOB_END_TIME --exe maia/src/maia

