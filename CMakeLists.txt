cmake_minimum_required(VERSION 3.14...3.16)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(SStar VERSION 0.1
              DESCRIPTION "Rewrite of SStar Estimator"
              LANGUAGES CXX)

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  message(STATUS "TESTING NOW")
  enable_testing()
  include(GoogleTest)
  add_subdirectory(tests)
endif()

macro(package_add_test TESTNAME FILES LIBRARIES)
    # create an exectuable in which the tests will be stored
    add_executable(${TESTNAME} ${FILES})
    # link the Google test infrastructure, mocking library,
    # and a default main fuction to the test executable.
    target_link_libraries(${TESTNAME} gtest gmock gtest_main ${LIBRARIES})
    # gtest_discover_tests replaces gtest_add_tests,
    gtest_discover_tests(${TESTNAME}
        WORKING_DIRECTORY ${PROJECT_DIR}
        PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "${PROJECT_DIR}"
        )
    set_target_properties(${TESTNAME} PROPERTIES FOLDER tests)
    target_compile_features(${TESTNAME} PUBLIC cxx_std_11)
    set_target_properties(${TESTNAME} PROPERTIES CXX_EXTENSIONS OFF)
    target_include_directories(${TESTNAME} PUBLIC include)
    add_test(NAME, ${TESTNAME} COMMAND ${TESTNAME})
endmacro()

# package_add_test(test1 test1.cpp lib_to_test)

# generate a lib
# add_library(lib STATIC src/lib.cpp include/sstar/lib.cpp)
# target_include_directoryies(lib PUBLIC include)
# target_compile_features(lib PUBLIC cxx_std_11)

# link to exe
# add_executable(target main.cpp)
# target_link_libraries(target PUBLIC lib)

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
  message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                 EXCLUDE_FROM_ALL)

add_library(vcf_file src/vcf_file.cpp include/sstar2/vcf_file.hpp)
target_include_directories(vcf_file PUBLIC include)
target_compile_features(vcf_file PUBLIC cxx_std_11)
set_target_properties(vcf_file PROPERTIES CXX_EXTENSIONS OFF)
package_add_test(vcf_file_test tests/test_vcf_file.cpp vcf_file)

# add_library(population_data src/population_data.cpp include/sstar2/population_data.hpp)
# target_include_directories(population_data PUBLIC include)
# target_compile_features(population_data PUBLIC cxx_std_11)
# package_add_test(population_data_test tests/test_population_data.cpp population_data)
