cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(RcppPlanc VERSION 1.0.0 LANGUAGES C CXX)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
include(CheckIncludeFile)
include(CheckCXXCompilerFlag)
include(cmake/CPM.cmake)
include(cmake/FindRModule.cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_C_FLAGS_RELEASE "")
set(CMAKE_CXX_FLAGS_RELEASE "")
set(CMAKE_VERBOSE_MAKEFILE TRUE)

if(APPLE)
    set(CMAKE_PROGRAM_PATH /opt/R/${CMAKE_HOST_SYSTEM_PROCESSOR}/bin)
    set(CMAKE_LIBRARY_PATH /opt/R/${CMAKE_HOST_SYSTEM_PROCESSOR}/lib)
    set(CMAKE_INCLUDE_PATH /opt/R/${CMAKE_HOST_SYSTEM_PROCESSOR}/include)
endif()

add_library(RcppPlanc SHARED "${PROJECT_SOURCE_DIR}/src/rcppplanc_nmf.cpp" "${PROJECT_SOURCE_DIR}/src/RcppExports.cpp" "${PROJECT_SOURCE_DIR}/src/conversion.cpp")
set_target_properties(RcppPlanc PROPERTIES PREFIX "")
if(APPLE)
set_target_properties(RcppPlanc PROPERTIES SUFFIX ".so")
endif()
target_include_directories(RcppPlanc
PRIVATE "${PROJECT_SOURCE_DIR}/src/common;${PROJECT_SOURCE_DIR}/src/nnls;${PROJECT_SOURCE_DIR}/src/nmf")

set(MKL_DIR "${PROJECT_SOURCE_DIR}/cmake")
find_package(MKL QUIET)
if(MKL_FOUND AND NOT DEFINED BLAS_VENDOR AND NOT WIN32)
    message(STATUS "found MKL")
    set(USE_MKL)
    set(BLAS_LIBRARIES "${MKL_LIBRARIES}")
    set(LAPACK_LIBRARIES "${MKL_LAPACK_LIBRARIES}")
    target_compile_definitions(RcppPlanc PRIVATE "MKL_FOUND")
    set(CBLAS_H_DIR "${MKL_INCLUDE_DIR}")
elseif(NOT WIN32)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif()

if(NOT(APPLE AND(CMAKE_HOST_SYSTEM_PROCESSOR MATCHES arm64)))
    find_package(OpenMP)
endif()

find_package(R REQUIRED)

find_r_module(Rcpp)
find_r_module(RcppArmadillo)
find_r_module(RcppProgress)

if(R_INCLUDE_DIR)
target_include_directories(RcppPlanc PRIVATE "${R_INCLUDE_DIR}")
else()
message(FATAL_ERROR "missing R headers")
endif()
if(R_RCPP)
target_include_directories(RcppPlanc PRIVATE "${R_RCPP}/include")
else()
message(FATAL_ERROR "missing RcppArmadillo")
endif()
if(R_RCPPARMADILLO)
target_include_directories(RcppPlanc PRIVATE "${R_RCPPARMADILLO}/include")
else()
message(FATAL_ERROR "missing RcppArmadillo")
endif()
if(R_RCPPPROGRESS)
target_include_directories(RcppPlanc PRIVATE "${R_RCPPPROGRESS}/include")
else()
message(FATAL_ERROR "missing RcppProgress")
endif()

add_library(hw_detect OBJECT "${PROJECT_SOURCE_DIR}/src/common/hw_detect.cpp")

if(WIN32)
    CPMAddPackage(NAME hwloc
        URL https://download.open-mpi.org/release/hwloc/v2.9/hwloc-2.9.3.tar.gz
        URL_HASH SHA256=5985db3a30bbe51234c2cd26ebe4ae9b4c3352ab788b1a464c40c0483bf4de59
        VERSION 2.9.3
        SOURCE_SUBDIR contrib/windows-cmake
        OPTIONS "HWLOC_ENABLE_TESTING OFF" "HWLOC_SKIP_LSTOPO ON" "HWLOC_SKIP_TOOLS ON"
    )
else()
    find_package(HWLOC)
    if(HWLOC_FOUND)
    target_link_libraries(hw_detect PRIVATE "${HWLOC_LDFLAGS}")
    target_include_directories(hw_detect PRIVATE "${HWLOC_INCLUDE_DIRS}")
    else()
        CPMAddPackage(NAME hwloc
            URL https://download.open-mpi.org/release/hwloc/v2.9/hwloc-2.9.3.tar.gz
            URL_HASH SHA256=5985db3a30bbe51234c2cd26ebe4ae9b4c3352ab788b1a464c40c0483bf4de59
            VERSION 2.9.3
            DOWNLOAD_ONLY
        )
    #find_package(Autotools REQUIRED)
    add_custom_command(OUTPUT "${hwloc_SOURCE_DIR}/Makefile"
    COMMAND ${hwloc_SOURCE_DIR}/configure --with-pic --prefix=${hwloc_SOURCE_DIR}/out/ --disable-picky --enable-shared=no --enable-static=yes
    --disable-cairo --disable-libxml2 --disable-io --disable-plugin-dlopen --disable-plugin-ltdl --enable-embedded-mode
    WORKING_DIRECTORY "${hwloc_SOURCE_DIR}"
    VERBATIM)
    add_custom_command(OUTPUT "${hwloc_SOURCE_DIR}/hwloc/.libs/libhwloc_embedded.a"
    DEPENDS "${hwloc_SOURCE_DIR}/Makefile"
    COMMAND "${CMAKE_MAKE_PROGRAM}"
    WORKING_DIRECTORY "${hwloc_SOURCE_DIR}")
    add_custom_target(hwloc DEPENDS "${hwloc_SOURCE_DIR}/hwloc/.libs/libhwloc_embedded.a")
    target_include_directories(hw_detect PRIVATE "${hwloc_SOURCE_DIR}/include")
    add_dependencies(hw_detect hwloc)
    target_link_libraries(hw_detect PRIVATE "${hwloc_SOURCE_DIR}/hwloc/.libs/libhwloc_embedded.a")
    if(APPLE)
    target_link_libraries(hwloc INTERFACE "-framework Foundation" "-framework IOKit")
    endif()
    endif()
endif()
if(WIN32)
    CPMAddPackage(NAME HighFive
        GIT_TAG v2.7.1
        GITHUB_REPOSITORY BlueBrain/HighFive
        OPTIONS "HIGHFIVE_EXAMPLES OFF" "HIGHFIVE_BUILD_DOCS OFF" "HIGHFIVE_USE_BOOST OFF"
        DOWNLOAD_ONLY TRUE)
    target_include_directories(RcppPlanc PRIVATE "${HighFive_SOURCE_DIR}/include")
    target_link_libraries(hw_detect PUBLIC hwloc)
    target_link_libraries(RcppPlanc PRIVATE -lhdf5 -lsz -lz -lopenblas -lgfortran -lgomp -lpthread -ldl -lm -lquadmath "${R_HOME}/bin/x64/R.dll")
else()
    CPMAddPackage(NAME HighFive
        GIT_TAG v2.7.1
        GITHUB_REPOSITORY BlueBrain/HighFive
        OPTIONS "HIGHFIVE_EXAMPLES OFF" "HIGHFIVE_BUILD_DOCS OFF" "HIGHFIVE_USE_BOOST OFF")
    target_link_libraries(RcppPlanc PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} HighFive)
endif()

target_link_libraries(RcppPlanc PRIVATE hw_detect)

if(OpenMP_CXX_FOUND)
target_compile_options(RcppPlanc PRIVATE ${OpenMP_CXX_FLAGS})
target_include_directories(RcppPlanc PRIVATE ${OpenMP_CXX_INCLUDE_DIRS})
endif()

check_cxx_compiler_flag("-funsafe-math-optimizations" UNSAFE_MATH)
if(UNSAFE_MATH)
    target_compile_options(RcppPlanc PRIVATE "-funsafe-math-optimizations")
endif()
check_cxx_compiler_flag("-ffinite-math-only" FINITE_MATH)
if(FINITE_MATH)
    target_compile_options(RcppPlanc PRIVATE "-ffinite-math-only")
endif()
check_cxx_compiler_flag("-fno-trapping-math" TRAPPING_MATH)
if(TRAPPING_MATH)
    target_compile_options(RcppPlanc PRIVATE "-fno-trapping-math")
endif()
check_cxx_compiler_flag("-fno-math-errno" ERRNO_MATH)
if(ERRNO_MATH)
    target_compile_options(RcppPlanc PRIVATE "-fno-math-errno")
endif()
check_cxx_compiler_flag("-fno-signed-zeros" UNSIGNED_ZERO)
if(UNSIGNED_ZERO)
    target_compile_options(RcppPlanc PRIVATE "-fno-signed-zeros")
endif()
check_cxx_compiler_flag("-freciprocal-math" RECIPROCAL_MATH)
if(RECIPROCAL_MATH)
    target_compile_options(RcppPlanc PRIVATE "-freciprocal-math")
endif()

set(CBLAS_H_DIR CACHE PATH "Directory containing cblas.h")

if(NOT APPLE AND NOT USE_MKL)
if(CBLAS_H_DIR)
target_include_directories(RcppPlanc PRIVATE ${CBLAS_H_DIR})
else()
check_include_file("cblas.h" CBLAS_HEADER)
if(NOT CBLAS_HEADER)
check_include_file("flexiblas/cblas.h" FLEXIBLAS_HEADER)
if(FLEXIBLAS_HEADER)
target_compile_definitions(RcppPlanc PRIVATE "HAVE_FLEXIBLAS_CBLAS_H")
else()
check_include_file("openblas/cblas.h" OPENBLAS_HEADER)
if(OPENBLAS_HEADER)
target_compile_definitions(RcppPlanc PRIVATE "HAVE_OPENBLAS_CBLAS_H")
else()
message(FATAL_ERROR "no cblas header, please set include directory with CBLAS_H_DIR variable")
endif()
endif()
endif()
endif()
endif()
