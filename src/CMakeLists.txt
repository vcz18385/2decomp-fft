if (BUILD_TARGET MATCHES "mpi_cpu")
  add_subdirectory(mpi_cpu)
endif (BUILD_TARGET MATCHES "mpi_cpu")

if (BUILD_TARGET MATCHES "mpi_gpu")
  add_subdirectory(mpi_gpu)
endif (BUILD_TARGET MATCHES "mpi_gpu")
