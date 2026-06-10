# Native (arm64 macOS, clang) build of the authors' tuttekiri C++ sources.
#
# The upstream tree baseline/tuttekiri/ is READ-ONLY: nothing here writes into
# it.  Its .cpp files are compiled in place; bind.cpp (the emscripten binding
# layer) is excluded and replaced by src/main.cpp + src/json_io.cpp.
cmake_minimum_required(VERSION 3.22)
project(tuttekiri_native CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POLICY_VERSION_MINIMUM 3.22)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(UPSTREAM ${CMAKE_CURRENT_SOURCE_DIR}/../tuttekiri/code/cpp)
if(NOT EXISTS ${UPSTREAM}/bind.cpp)
  message(FATAL_ERROR "upstream sources not found at ${UPSTREAM}")
endif()

include(FetchContent)

# --- Eigen ------------------------------------------------------------------
# Homebrew ships Eigen 5.0.1, which libigl v2.5.0 and Optiz do not support
# (see baseline/README.md).  Pin the same 3.4.0 the authors' CMakeLists falls
# back to, so the numerics match theirs bit for bit.
FetchContent_Declare(Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  SOURCE_SUBDIR cmake)          # same hack the authors use to skip Eigen's CMakeLists
FetchContent_MakeAvailable(Eigen3)
if(NOT TARGET Eigen3::Eigen)
  add_library(Eigen3::Eigen INTERFACE IMPORTED)
  set_target_properties(Eigen3::Eigen PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${eigen3_SOURCE_DIR})
endif()

# --- libigl + Optiz (exactly the authors' pins) ------------------------------
set(LIBIGL_RESTRICTED_TRIANGLE ON CACHE BOOL "" FORCE)
FetchContent_Declare(libigl
  GIT_REPOSITORY https://github.com/libigl/libigl.git
  GIT_TAG v2.5.0)
# Optiz's CMakeLists has an APPLE block that appends a bare `-fopenmp` whenever
# the platform is Apple, guarded by CMAKE_CXX_COMPILER_ID STREQUAL "Clang".
# Apple's compiler reports "AppleClang", so the `-Xclang` prefix is skipped and
# Apple clang rejects the bare flag.  Optiz's own find_package(OpenMP) has
# already failed on this machine (no libomp), so the intended configuration is
# no-OpenMP; the patch below just deletes those two lines.  Single-threaded
# Optiz is also bit-reproducible, which is what a referee baseline wants.
FetchContent_Declare(Optiz
  GIT_REPOSITORY https://github.com/segaviv/optiz.git
  GIT_TAG master
  PATCH_COMMAND sed -i "" -e "/-Xclang/d" -e "/-fopenmp/d" CMakeLists.txt
  UPDATE_DISCONNECTED TRUE)
FetchContent_MakeAvailable(Optiz libigl)

# --- upstream sources, minus bind.cpp ---------------------------------------
file(GLOB UPSTREAM_SRC
  ${UPSTREAM}/utils/*.cpp
  ${UPSTREAM}/coloring/*.cpp
  ${UPSTREAM}/geometry/*.cpp
  ${UPSTREAM}/parameterization/*.cpp
  ${UPSTREAM}/opt/*.cpp)

# --- WP5 crash fix (opt-in, -DKIRI_CRASHFIX=ON) ------------------------------
# Two upstream .cpp files are compiled from patched copies in src_fixed/ instead
# of in place (src_fixed/*.cpp = upstream + results/kill/native200/crashfix.patch;
# baseline/tuttekiri/ itself stays byte-identical).  Build it into a SEPARATE
# directory so build/tuttekiri_cli, the binary the 600 s run used, is untouched:
#   cmake -S baseline/native -B baseline/native/build_fixed -DKIRI_CRASHFIX=ON \
#         -DCMAKE_OSX_ARCHITECTURES=arm64 && cmake --build baseline/native/build_fixed
option(KIRI_CRASHFIX "compile the F24-patched copies in src_fixed/" OFF)
if(KIRI_CRASHFIX)
  foreach(f ${UPSTREAM_SRC})
    if(f MATCHES "utils/Hmesh\\.cpp$" OR f MATCHES "geometry/unit_pattern\\.cpp$")
      list(REMOVE_ITEM UPSTREAM_SRC ${f})
    endif()
  endforeach()
  list(APPEND UPSTREAM_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/src_fixed/Hmesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src_fixed/unit_pattern.cpp)
endif()

add_executable(tuttekiri_cli
  ${UPSTREAM_SRC}
  src/json_io.cpp
  src/their_system.cpp
  src/bind_port.cpp
  src/main.cpp)

# emsc_stub must come first so #include <emscripten/bind.h> resolves to it.
target_include_directories(tuttekiri_cli BEFORE PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/emsc_stub)
target_include_directories(tuttekiri_cli PRIVATE
  ${UPSTREAM}
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  /opt/homebrew/include)          # nlohmann/json
if(KIRI_CRASHFIX)
  # the relocated copies still include their neighbours by bare name
  target_include_directories(tuttekiri_cli PRIVATE
    ${UPSTREAM}/utils ${UPSTREAM}/geometry)
endif()

target_link_libraries(tuttekiri_cli PRIVATE
  Eigen3::Eigen Optiz igl::core igl_restricted::triangle)

target_compile_options(tuttekiri_cli PRIVATE
  -Wno-deprecated-declarations -Wno-return-type -Wno-unused-value)
