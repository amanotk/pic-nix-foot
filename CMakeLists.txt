cmake_minimum_required(VERSION 3.10)

# PIC-NIX source tree
if(NOT DEFINED PICNIX_DIR)
  message(FATAL_ERROR "Please define PICNIX_DIR")
endif()
add_compile_definitions(PICNIX_DIR="${PICNIX_DIR}")

# common configuration
include(${PICNIX_DIR}/common.cmake OPTIONAL)

# project
project(foot)

# source directory
add_subdirectory(${PICNIX_DIR}/nix ${CMAKE_BINARY_DIR}/nix)
add_subdirectory(${PICNIX_DIR}/expic ${CMAKE_BINARY_DIR}/expic)

# include directory
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${PICNIX_DIR})
include_directories(${PICNIX_DIR}/nix)
include_directories(${PICNIX_DIR}/nix/thirdparty)
include_directories(${PICNIX_DIR}/nix/expic)

# build executable
add_executable(foot.out main.cpp)

# link library
target_link_libraries(foot.out nix expic)
