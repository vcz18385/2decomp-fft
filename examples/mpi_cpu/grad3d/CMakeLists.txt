file(GLOB files_test grad3d.f90)
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(grad3d ${files_test})
target_link_libraries(grad3d PRIVATE decomp2d_cpu examples_utils)

# Run the test(s)
set(run_dir "${test_dir}/grad3d")
file (COPY "${CMAKE_SOURCE_DIR}/examples/mpi_cpu/grad3d/adios2_config.xml" DESTINATION ${run_dir})
message(STATUS "Example dir ${run_dir}")
file(MAKE_DIRECTORY ${run_dir})
add_test(NAME grad3d COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:grad3d> ${TEST_ARGUMENTS} WORKING_DIRECTORY ${run_dir})
