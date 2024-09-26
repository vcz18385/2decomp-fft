export FC="mpif90"

rm -rf `pwd`/build_mpi_cpu
cmake -S `pwd` -B `pwd`/build_mpi_cpu -DBUILD_TARGET=mpi_cpu -DFFT_Choice=generic -DBUILD_TESTING=ON
cd `pwd`/build_mpi_cpu; make; cd -

rm -rf `pwd`/build_mpi_gpu
cmake -S `pwd` -B `pwd`/build_mpi_gpu -DBUILD_TARGET=mpi_gpu -DFFT_Choice=cufft -DBUILD_TESTING=ON -DENABLE_OPENACC=ON
cd `pwd`/build_mpi_gpu; make; cd -

mkdir `pwd`/build
mkdir `pwd`/build/lib

cp `pwd`/build_mpi_cpu/lib/lib* `pwd`/build/lib/
cp `pwd`/build_mpi_gpu/lib/lib* `pwd`/build/lib/

