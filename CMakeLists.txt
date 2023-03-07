cmake_minimum_required(VERSION 3.10)

# PIC-NIX source tree in absolute path
if(NOT DEFINED PICNIX_DIR)
  message(FATAL_ERROR "Please define PICNIX_DIR")
endif()
get_filename_component(PICNIX_DIR ${PICNIX_DIR} ABSOLUTE)

# common configuration
include(${PICNIX_DIR}/common.cmake OPTIONAL)

set(PROJECT_NAME "foot")

# project
project(${PROJECT_NAME} CXX)

# source directory
add_subdirectory(${PICNIX_DIR}/nix ${CMAKE_BINARY_DIR}/nix)
add_subdirectory(${PICNIX_DIR}/expic ${CMAKE_BINARY_DIR}/expic)

# executable
add_executable(${PROJECT_NAME} main.cpp)
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME main.out)

# include directory
target_include_directories(${PROJECT_NAME} PRIVATE ${PICNIX_DIR})

# link library
target_link_libraries(${PROJECT_NAME} nix expic)
