rm -rf `pwd`/build_mpi_NVHPC241_openacc_cuda
export FC="mpif90"
cmake -S `pwd` -B `pwd`/build_mpi_NVHPC241_openacc_cuda -DBUILD_TARGET=mpi_gpu -DFFT_Choice=cufft -DBUILD_TESTING=ON -DENABLE_OPENACC=ON
