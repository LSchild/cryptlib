cmake_minimum_required(VERSION 3.26)
project(large_functions)

set(CMAKE_CXX_STANDARD 17)
include(FetchContent)

# SETUP OpenFHE

find_package(OpenFHE CONFIG REQUIRED PATHS "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/openfhe" NO_DEFAULT_PATH)
# find_package(OpenFHE CONFIG REQUIRED)

if (OpenFHE_FOUND)
    message(STATUS "FOUND PACKAGE OpenFHE")
    message(STATUS "OpenFHE Version: ${BASE_OPENFHE_VERSION}")
    message(STATUS "OpenFHE installed as shared libraries: ${OpenFHE_SHARED}")
    message(STATUS "OpenFHE include files location: ${OpenFHE_INCLUDE}")
    message(STATUS "OpenFHE lib files location: ${OpenFHE_LIBDIR}")
    message(STATUS "OpenFHE Native Backend size: ${OpenFHE_NATIVE_SIZE}")
else()
    message(FATAL_ERROR "PACKAGE OpenFHE NOT FOUND")
endif ()

#set( CMAKE_CXX_FLAGS "${OpenFHE_CXX_FLAGS} -flto" )
set(CMAKE_CXX_FLAGS "-march=native -g2 ")

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${OpenFHE_INCLUDE} )
include_directories( ${OpenFHE_INCLUDE}/third-party/include )
include_directories( ${OpenFHE_INCLUDE}/core )
include_directories( ${OpenFHE_INCLUDE}/pke )
include_directories( ${OpenFHE_INCLUDE}/binfhe )

### add directories for other OpenFHE modules as needed for your project

link_directories( ${OpenFHE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
if(BUILD_STATIC)
    set( CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
    link_libraries( ${OpenFHE_STATIC_LIBRARIES} )
else()
    set( CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS} )
    link_libraries( ${OpenFHE_SHARED_LIBRARIES} )
endif()

# SETUP HEXL
set(HEXL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/hexl")
include_directories("${HEXL_DIR}/include")

# SETUP GTEST
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# SETUP OMP
#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

add_library(large_functions SHARED
        lib/include/rgsw.h
        lib/src/rgsw.cpp
        lib/include/automorphism_key.h
        lib/src/automorphism_key.cpp
        lib/include/rlwe_key_switching.h
        lib/src/rlwe_key_switching.cpp
        lib/include/utils.h
        lib/src/utils.cpp
        lib/include/setup.h
        lib/include/blind_rotator.h
        lib/src/blind_rotator.cpp
        lib/include/scheme_switch.h
        lib/src/scheme_switch.cpp
        lib/include/digit_decomposer.h
        lib/src/digit_decomposer.cpp
        lib/include/functional_bootstrap.h
        lib/src/functional_bootstrap.cpp
        lib/include/giant_lut_evaluator.h
        lib/src/giant_lut_evaluator.cpp
        lib/include/lwe_key_switching.h
        lib/src/lwe_key_switching.cpp
        lib/include/crt_database.h
        lib/src/crt_database.cpp
        lib/include/response_compression.h
        lib/src/response_compression.cpp
        lib/include/pirouette.h
        lib/src/pirouette.cpp
        lib/src/default_params.cpp
        lib/include/digit_decomposer_fast.h
        lib/src/digit_decomposer_fast.cpp
        lib/include/interfaces.h
        lib/include/cggi_blind_rotator.h
        lib/include/bmmp_blind_rotator.h
        lib/include/mux_operator.h
        lib/src/mux_operator.cpp
        lib/src/cggi_blind_rotator.cpp
        lib/include/base_crypto.h
        lib/src/base_crypto.cpp
        lib/src/bmmp_blind_rotator.cpp
        lib/include/speed_utils.h
        lib/include/glwe_conversion.h
        lib/src/lwe_conversion.cpp
        lib/include/gadget_decomp.h
        lib/src/rlwe_conversion.cpp
        lib/include/math_utils.h
        lib/src/math_utils.cpp
        lib/include/rgsw_conversion.h
        lib/include/container_types.h
)

target_include_directories(large_functions PUBLIC lib/include)
target_link_libraries(large_functions PUBLIC ${HEXL_DIR}/lib64/libhexl.a)

# Testing
enable_testing()

add_executable(test_library tests/dummy_test.cpp
        tests/rlwe_tests.cpp
        tests/scheme_switch_tests.cpp
        tests/blind_rotator_tests.cpp
        tests/digit_decomp_test.cpp
        tests/lwe_keyswitch_tests.cpp
        tests/Giant_lut_eval_test.cpp
        tests/pirouette_test.cpp
        tests/digit_decomp_tests_fast.cpp
        tests/mux_operator_tests.cpp
        tests/base_crypto_tests.cpp
        tests/cggi_blind_rotation_test.cpp
        tests/bmmp_blind_rotation_tests.cpp
        tests/glwe_conversion_test.cpp
        tests/noise_tests.cpp
        tests/automorphism_test.cpp)
target_link_libraries(test_library
        large_functions
    GTest::gtest_main
)
target_include_directories(test_library PUBLIC lib/include)

include(GoogleTest)
gtest_discover_tests(test_library)

add_executable(run_pirouette run_pirouette.cpp)
target_include_directories(run_pirouette PUBLIC lib/include)
target_link_libraries(run_pirouette large_functions)

#add_executable(pirouette_encrypt iot/pirouette_encrypt.cpp)
#add_executable(respire_encrypt iot/respire_transciphering_encrypt.cpp)