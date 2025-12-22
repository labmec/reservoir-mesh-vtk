cmake_minimum_required(VERSION 3.14)
project(TestErick)

# Prefer a user-provided NeoPZ_DIR or add common install locations to CMAKE_PREFIX_PATH
if(NOT DEFINED NeoPZ_DIR)
  set(NeoPZ_DIR "${CMAKE_SOURCE_DIR}/../../../firstPZ/neopz_install" CACHE PATH "NeoPZ installation directory")
endif()
list(APPEND CMAKE_PREFIX_PATH "${NeoPZ_DIR}")

find_package(NeoPZ REQUIRED)

add_custom_target( GlocalFiles Sources ${PZ_SRC})

add_executable(TestErick meshVisualization.cpp)

target_link_libraries(TestErick NeoPZ::pz)