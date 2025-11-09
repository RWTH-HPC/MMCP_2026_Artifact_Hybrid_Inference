Point data feature in the DG block
==================================

Introduction
------------

The point data feature records the state of the conservative variables at
specified intervals at specified points. The points can be anywhere in the
domain you're computing in. Some restrictions apply to the sample interval.


Format
------

All (input and output) text files are in the gnuplot standard format. That
means one data point per line and the values are separated by spaces.


Usage
-----

To use the point data feature you first have to create an input file which
contains the points you are interested in. The points can have any shape you
want.
If you just want a straight line you can use the script under
`preprocessing/gen_line_points_csv.py` to generate them (Use -h option to the
help).
For a circular shape see the script under `preprocessing/gen_circle_points.py`.
After that you specify the property `pointDataFileName_0` (which is a string)
with the path to your input file.
In addition you also need to specify two intervals:
The first one is the sample interval (`pointDataSampleInterval`) which
specifies how often the state of the variables should be recorded. The default
value is 1.
The second one is the write interval (`pointDataWriteInterval`) which specifies
how often the recorded values should be written out to a file.
The point of the write interval is to minimize expensive I/O time. So choose a
large enough value. Be aware that all the data which is sampled but not written
to disk is of course buffered in RAM. So don't overdo it. Also be aware of the
restrictions explained below.
The write interval defaults to the restart interval, so you don't need to
specify it, if you specify a restart interval. To be able to restart a
simulation without data loss the write interval needs to be a divisor of the
restart interval.
The sample interval also needs to be a divisor of the write interval for
technical reasons. Also it makes much more sense this way.
You can also set the start and end time step of recording by the properties
`pointDataStartTimeStep` and `pointDataEndTimeStep`. In default case the
recording will start at time step 0 and continues to the last timestep of
the simulation.
To process multiple input files just add another `pointDataFileName_x` to the
properties, where x starts at 0 and increases with 1 for each additional file.
You can set a property for all input files when you just set a property without
the number extension. To specify one input file, add a property with a number
extension at the end related to the input file. 
Finally you can optionally declare `pointDataEnabled` to have an easy way to
switch the entire feature on and off.

After MAIA is finished you should have several (at least one)
`point_data_*.Netcdf` files in your `out` folder.

Since the I/O is optimized for performance the files are kind of hard to use.
Therefore it is recommended to use the
`postprocessing/dg_process_point_data.py` script to get the recorded data in a
meaningful format.


Full Example
------------

Assume a 2D testcase in which you are interested in the states along a line.
You want to sample every 5 timesteps. The line goes from (-20, 1e-10) to (20,
1e-10). You want to have 100 points along this line.

Generate line:
```
~/code/trunk/preprocessing/gen_line_points.py 2 100 -20 1e-10 20 1e-10 -o line
```

Add the following properties to your properties_run.cdl file (Don't forget to
edit the header):
```
pointDataFileName_0 = "line"
pointDataSampleInterval = 5
pointDataWriteInterval = 100
```

Run MAIA.
Extract data for points 5 to 20 from time 0.7 to 2.1:
```
~/code/trunk/postprocessing/process_point_data.py -p 5:20 -t 0.7:2.1 -i .
```

Add the following properties to your properties_run.cdl file to process
e.g. two input files with different settings.
```
pointDataFileName_0 = "line";
pointDataFileName_1 = "circle";
pointDataSampleInterval = 5;
pointDataWriteInterval_0 = 100;
pointDataWriteInterval_1 = 50;
```

This will create two output file series, where both where sampled every fifth
time step, but there will be twice as many circle files than line files.
