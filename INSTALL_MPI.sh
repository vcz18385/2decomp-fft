rm -rf `pwd`/build_mpi_cpu
export FC="mpif90"
cmake -S `pwd` -B `pwd`/build_mpi_cpu -DBUILD_TARGET=mpi_cpu -DFFT_Choice=generic -DBUILD_TESTING=ON
