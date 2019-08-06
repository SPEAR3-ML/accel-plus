cmake_minimum_required(VERSION 3.15.0)
project(accel+)
# The version number.
set(accel+_VERSION_MAJOR 1)
set(accel+_VERSION_MINOR 0)

# Disable in-source builds to prevent source tree corruption.
if("${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "
FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.
")
endif()

# Build shared libs or static libs
option(BUILD_SHARED_LIBS "Build shared libs" OFF)
option(WINDOWS_EXPORT_ALL_SYMBOLS "Auto export all symbols in shared libs" ON)

# Setup the lib and bin directory
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

include_directories("${PROJECT_SOURCE_DIR}/lib/accel")
add_subdirectory(lib/accel)

set(SRCS process.cpp)

add_executable(accel+ ${SRCS})
target_link_libraries(accel+ accel)
