#!/usr/local_rwth/bin/zsh

source ./setup_env_claix23.sh

source ./cpp-ml-interface/extern/python/venv/bin/activate

export SCOREP_WRAPPER_INSTRUMENTER_FLAGS="--verbose=1 --nocompiler --user --mpp=mpi --io=none --memory=none --thread=none --nocuda"
export SCOREP_ENABLE_CUDA=0
#export SCOREP_WRAPPER_COMPILER_FLAGS="-g -DSCOREP"
pushd maia

./configure.py 1 2 --enable-instrumentation scorep --instrument mpi --instrument user && make -j