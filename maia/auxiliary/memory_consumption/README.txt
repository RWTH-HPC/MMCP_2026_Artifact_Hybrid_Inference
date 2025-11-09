This program allocates a chunk of memory and then displays the total memory allocated (Virtual memory)
and the current memory usage (RAM).

To compile and run the program on the AIA cluster:
$ mpic++ memoryConsumption.cpp
$ mpirun -np <number_of_processes> ./a.out <memory size per process in MB>
