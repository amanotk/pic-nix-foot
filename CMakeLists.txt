cmake_minimum_required(VERSION 3.14)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE TRUE)

# project
project(pic-nix-foot CXX)

# pic-nix as an external project
include(FetchContent)

set(PICNIX_DIR ${CMAKE_CURRENT_BINARY_DIR}/pic-nix)
FetchContent_Declare(
  pic-nix
  GIT_REPOSITORY https://github.com/amanotk/pic-nix.git
  GIT_TAG main
  SOURCE_DIR ${PICNIX_DIR}
  BINARY_DIR ${PICNIX_DIR}-build
)
set(BUILD_EXAMPLE OFF) # disable building example
FetchContent_MakeAvailable(pic-nix)

# Filesystem package
set(CMAKE_MODULE_PATH ${pic-nix_SOURCE_DIR}/nix)
find_package(Filesystem REQUIRED)

# order of shape function
set(SHAPE_ORDER 2 CACHE STRING "Order of shape function")
add_compile_definitions(SHAPE_ORDER=${SHAPE_ORDER})

# executable
add_executable(${PROJECT_NAME} main.cpp)
set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME main.out)

# include directory
target_include_directories(${PROJECT_NAME} PRIVATE ${PICNIX_DIR})

# link library
target_link_libraries(${PROJECT_NAME} nix)
target_link_libraries(${PROJECT_NAME} expic)
target_link_libraries(${PROJECT_NAME} std::filesystem)

# copy
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/alfven DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/buneman DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
