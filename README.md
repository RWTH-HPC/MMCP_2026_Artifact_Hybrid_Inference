# Hybrid Inference Optimization for AI-Enhanced Turbulent Boundary Layer Simulation on Heterogeneous Systems - Supplemental Material

This is supplemental material for the paper "Hybrid Inference Optimization for AI-Enhanced Turbulent Boundary Layer Simulation on Heterogeneous Systems" submitted to the MMCP 2026 workshop.

!!! Please note that this repository will be updated until the conference date January 26-29, 2026 !!!

The code was tested on the CLAIX-2023 cluster at RWTH Aachen University using the following software stack
```
module load foss/2023b
module load FFTW.MPI/3.3.10
module load PnetCDF/1.12.3
module load HDF5/1.14.3
module load Python/3.11.5
module load CMake/3.29.3
module load imkl/2024.2.0

module load CUDA/12.4.0
module load cuDNN/8.9.7.29-CUDA-12.4.0
module load Score-P/8.4-CUDA-12.4.0
```
Please report any issue on other systems!

## Repository Structure
- [cpp-ml-interface](cpp-ml-interface/): Source code of the ML-Module partially ported to C++ from [Fortran](https://github.com/RWTH-HPC/Fortran-ML-Interface). It also includes the [AIxeleratorService](https://github.com/RWTH-HPC/AIxeleratorService/tree/MMCP_2026) library for hybrid inference as a submodule.
- [maia](maia/): Source code of the finite volume solver from the multi-physics simulation framework [m-AIA](https://git.rwth-aachen.de/aia/m-AIA/m-AIA/-/tree/v2024.1?ref_type=tags) coupled with the TBL-Transformer model.
- [input](input/): contains the requird input files for the solver (e.g. grid, restart file, Transformer model)