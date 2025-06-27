cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 17)

# project
set(PROJECT_NAME "foot")
project(${PROJECT_NAME} CXX)

# xsimd for xtensor
add_compile_definitions(XTENSOR_USE_XSIMD)

# MPI_THREAD_MULTIPLE
option(MPI_THREAD_MULTIPLE "Enable MPI_THREAD_MULTIPLE" ON)
if(MPI_THREAD_MULTIPLE)
    add_compile_definitions(HAS_MPI_THREAD_MULTIPLE)
endif()

# pic-nix as an external project
include(FetchContent)

set(PICNIX_BUILD_EXAMPLE OFF) # disable building example

# pic-nix directory
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/pic-nix")
  # use local pic-nix directory if exists
  set(PICNIX_DIR "${CMAKE_CURRENT_SOURCE_DIR}/pic-nix")
  add_subdirectory(${PICNIX_DIR} ${CMAKE_BINARY_DIR}/pic-nix)
  message(STATUS "Using local pic-nix directory : ${PICNIX_DIR}")
else()
  # download pic-nix from github
  set(PICNIX_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pic-nix)
  FetchContent_Declare(
    pic-nix
    GIT_REPOSITORY https://github.com/amanotk/pic-nix.git
    GIT_TAG v0.1.3
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/pic-nix
    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/pic-nix
  )
  FetchContent_MakeAvailable(pic-nix)
  message(STATUS "Using pic-nix directory cloned from github")
endif()

# Filesystem package
set(CMAKE_MODULE_PATH ${PICNIX_DIR}/nix)
find_package(Filesystem REQUIRED)

# executable
add_executable(${PROJECT_NAME} main.cpp)
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME main.out)

# include directory
target_include_directories(${PROJECT_NAME} PRIVATE ${PICNIX_DIR})

# link library
target_link_libraries(${PROJECT_NAME} std::filesystem)
target_link_libraries(${PROJECT_NAME} pic nix)

# copy
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/alfven DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/buneman DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
