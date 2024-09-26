file(GLOB files_test init_test.f90)
include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable(init_test ${files_test})
target_link_libraries(init_test PRIVATE decomp2d_cpu examples_utils)

# Run the test(s)
set(run_dir "${test_dir}/init_test")
message(STATUS "Example dir ${run_dir}")
file(MAKE_DIRECTORY ${run_dir})
add_test(NAME init_test COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} $<TARGET_FILE:init_test> ${TEST_ARGUMENTS} WORKING_DIRECTORY ${run_dir})
