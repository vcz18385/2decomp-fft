file(GLOB files_decomp factor.f90
	               decomp_2d_constants.f90
		       decomp_2d_mpi.f90
	               decomp_2d.f90
		       log.f90
                       ../io_utilities.f90
                       ../io_object_mpi.f90
                       ../io_object_adios.f90
                       ../io_mpi.f90
                       ../io_adios_none.f90
                       ../io_adios.f90
		       ../interp.f90
		       transpose_x_to_y.f90
		       transpose_y_to_x.f90
		       transpose_y_to_z.f90
		       transpose_z_to_y.f90)
include_directories(${CMAKE_SOURCE_DIR}/src)

if(FFTW_FOUND)
  if (${FFT_Choice} MATCHES "fftw_f03")
    message (STATUS "Compiling using FFTW3 F2003 interface")
    file(GLOB files_fft fft_fftw3_f03.f90)
  else()
    message (STATUS "Compiling using FFTW3 legacy interface")
    file(GLOB files_fft fft_fftw3.f90)
  endif()
elseif(MKL_FOUND)
  message (STATUS "Compiling using MKL")
  file(GLOB files_fft $ENV{MKLROOT}/include/mkl_dfti.f90 fft_mkl.f90)
elseif(CUFFT_FOUND)
  message (STATUS "Compiling using cuFFT")
  file(GLOB files_fft fft_cufft.f90)
else(FFTW_FOUND)
  message (STATUS "Compiling using Generic FFT")
  file(GLOB files_fft glassman.f90 fft_generic.f90)
endif(FFTW_FOUND)
list(APPEND files_fft fft_log.f90)

if (ENABLE_PROFILER)
  file(GLOB prof_files ../profiler_${ENABLE_PROFILER}.f90)
else (ENABLE_PROFILER)
  file(GLOB prof_files ../profiler_none.f90)
endif()

if (IO_BACKEND MATCHES "adios2")
  file(GLOB files_adios ../io_adios.f90 io_object_adios.f90)
else (IO_BACKEND MATCHES "adios2")
  file(GLOB files_adios ../io_adios_none.f90)
endif (IO_BACKEND MATCHES "adios2")

set(SRCFILES ${files_decomp} ${files_fft} ${prof_files} ${files_adios})

add_library(decomp2d_cpu ${SRCFILES})
target_include_directories(decomp2d_cpu PUBLIC
  $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
if (MPI_FOUND)
  target_link_libraries(decomp2d_cpu PRIVATE MPI::MPI_Fortran)
endif (MPI_FOUND)

if(FFTW_FOUND)
  if (FFTW_DOUBLE_LIB_FOUND)
    message(STATUS "FFTW_DOUBLE     : ${FFTW_DOUBLE_LIB}")
    target_link_libraries(decomp2d_cpu PRIVATE ${FFTW_DOUBLE_LIB})
  endif (FFTW_DOUBLE_LIB_FOUND)
  if (FFTW_FLOAT_LIB_FOUND)
    message(STATUS "FFTW_SINGLE     : ${FFTW_FLOAT_LIB}")
    target_link_libraries(decomp2d_cpu PRIVATE ${FFTW_FLOAT_LIB})
  endif (FFTW_FLOAT_LIB_FOUND)
  if (NOT (FFTW_DOUBLE_LIB_FOUND OR FFTW_FLOAT_LIB_FOUND))
    message(FATAL_ERROR "CMake step for FFTW failed")
  endif (NOT (FFTW_DOUBLE_LIB_FOUND OR FFTW_FLOAT_LIB_FOUND))
  target_include_directories(decomp2d_cpu PRIVATE ${FFTW_INCLUDE_DIRS})
elseif(MKL_FOUND)
  target_compile_options(decomp2d_cpu PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
  target_include_directories(decomp2d_cpu PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  target_link_libraries(decomp2d_cpu PUBLIC $<LINK_ONLY:MKL::MKL>)
endif(FFTW_FOUND)

if (ENABLE_PROFILER)
  target_compile_definitions(decomp2d_cpu PUBLIC -DPROFILER)
  if (caliper_FOUND)
    target_include_directories(decomp2d_cpu PUBLIC ${caliper_INCLUDE_DIR}/caliper/fortran)
    target_link_libraries(decomp2d_cpu PRIVATE caliper)
  endif()
endif()

if (HALO_DEBUG)
  target_compile_definitions(decomp2d_cpu PUBLIC -DHALO_DEBUG)
endif()

if (IO_BACKEND MATCHES "adios2")
  target_link_libraries(decomp2d_cpu PUBLIC adios2::fortran_mpi adios2::fortran)
endif (IO_BACKEND MATCHES "adios2")

install(TARGETS decomp2d_cpu
  EXPORT decomp2d_cpu-targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.mod")
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.smod")
      
## Packaging
install(EXPORT decomp2d_cpu-targets
  DESTINATION ${CMAKE_INSTALL_LIBDIR})

configure_file(
  ${CMAKE_SOURCE_DIR}/cmake/decomp2d_cpu-config.cmake.in
  ${PROJECT_BINARY_DIR}/cmake/decomp2d_cpu-config.cmake @ONLY
)
configure_file(
  ${CMAKE_SOURCE_DIR}/cmake/decomp2d_cpu-config-version.cmake.in
  ${PROJECT_BINARY_DIR}/cmake/decomp2d_cpu-config-version.cmake
)

install(FILES
  ${PROJECT_BINARY_DIR}/cmake/decomp2d_cpu-config.cmake
  ${PROJECT_BINARY_DIR}/cmake/decomp2d_cpu-config-version.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/decomp2d_cpu)

# Make 2decomp usable from build tree.
export(TARGETS decomp2d_cpu
  FILE ${PROJECT_BINARY_DIR}/decomp2d_cpu-targets.cmake)
configure_file(
  ${CMAKE_SOURCE_DIR}/cmake/decomp2d_cpu-build_config.cmake.in
  ${PROJECT_BINARY_DIR}/decomp2d_cpu-config.cmake @ONLY
)
