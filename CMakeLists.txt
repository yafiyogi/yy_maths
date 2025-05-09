#
#
#  MIT License
#
#  Copyright (c) 2025 Yafiyogi
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.
#
#

# PKG_CONFIG_PATH=${HOME}/code/tp/lib/pkgconfig:${HOME}/code/tp/share/pkgconfig:${PKG_CONFIG_PATH} cmake --fresh -DCMAKE_INSTALL_PREFIX=${HOME}/code/cpp/build -DYY_THIRD_PARTY_LIBRARY=${HOME}/code/tp -DCMAKE_POSITION_INDEPENDENT_CODE=TRUE -DCMAKE_BUILD_TYPE=Release .
# cmake --fresh -DCMAKE_INSTALL_PREFIX=${HOME}/code/cpp/build -DYY_THIRD_PARTY_LIBRARY=${HOME}/code/tp -DCMAKE_POSITION_INDEPENDENT_CODE=TRUE -DCMAKE_BUILD_TYPE=Release ..

cmake_minimum_required(VERSION 3.25)

project(yy_maths_proj LANGUAGES CXX VERSION 0.0.1)

set(CMAKE_EXPORT_COMPILE_COMMANDS)
set(YY_ROOT "${PROJECT_SOURCE_DIR}/..")
set(YY_CMAKE "${YY_ROOT}/yy_cmake")

include(${YY_CMAKE}/cmake_common.txt)

include(${YY_CMAKE}/cmake_boost-1.84.0.txt)
include(${YY_CMAKE}/cmake_clang_tidy.txt)

#find_package(fmt REQUIRED)

add_library(yy_maths STATIC)

target_compile_options(yy_maths
  PRIVATE
  "-DSPDLOG_COMPILED_LIB"
  "-DSPDLOG_FMT_EXTERNAL")

target_include_directories(yy_maths
  SYSTEM PRIVATE
    "${YY_THIRD_PARTY_LIBRARY}/include" )

target_sources(yy_maths
  PRIVATE
    yy_ekf.cpp
  PUBLIC FILE_SET HEADERS
    FILES
      yy_ekf.hpp
      yy_fib.hpp
      yy_diagonal_matrix.hpp
      yy_matrix.hpp
      yy_matrix_fmt.hpp
      yy_matrix_fwd.hpp
      yy_matrix_util.hpp)

install(TARGETS yy_maths
  EXPORT yy_mathsTargets
  FILE_SET HEADERS DESTINATION include/yy_maths)

install(EXPORT yy_mathsTargets
  NAMESPACE yy_maths::
  DESTINATION lib/cmake/yy_maths)

include(CMakePackageConfigHelpers)
write_basic_package_version_file("yy_mathsConfigVersion.cmake"
  VERSION ${yy_maths_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(FILES "yy_mathsConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/yy_mathsConfigVersion.cmake"
  DESTINATION lib/cmake/yy_maths)

#add_subdirectory(unit_tests)
add_subdirectory(examples)
#add_subdirectory(benchmarks)

add_yy_tidy_all(yy_maths)
