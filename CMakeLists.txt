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

# order of shape function
option(SHAPE_ORDER "Set SHAPE_ORDER to specify order of shape function" 1)
add_compile_definitions(SHAPE_ORDER=${SHAPE_ORDER})

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

# copy
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/alfven DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
