cmake_minimum_required(VERSION 3.14)
project(TestErick)

# Prefer a user-provided NeoPZ_DIR or add common install locations to CMAKE_PREFIX_PATH
if(NOT DEFINED NeoPZ_DIR)
  set(NeoPZ_DIR "${CMAKE_SOURCE_DIR}/../../../firstPZ/neopz_install" CACHE PATH "NeoPZ installation directory")
endif()
list(APPEND CMAKE_PREFIX_PATH "${NeoPZ_DIR}")

find_package(NeoPZ REQUIRED)

add_custom_target( GlocalFiles Sources ${PZ_SRC})

# # ------------------Executables------------------
# # add_executable(TestErick
# #     meshVisualization.cpp
# #     DarcySpaces.cpp
# #     mainDarcy.cpp
# # )
# add_executable(TestErick mainDarcy.cpp)

cmake_minimum_required(VERSION 3.14)
project(TestErick)

# Prefer a user-provided NeoPZ_DIR or add common install locations to CMAKE_PREFIX_PATH
if(NOT DEFINED NeoPZ_DIR)
  set(NeoPZ_DIR "${CMAKE_SOURCE_DIR}/../../../firstPZ/neopz_install" CACHE PATH "NeoPZ installation directory")
endif()
list(APPEND CMAKE_PREFIX_PATH "${NeoPZ_DIR}")

find_package(NeoPZ REQUIRED)

# add_custom_target( GlocalFiles Sources ${PZ_SRC})


# # ------------------Executables------------------
# # add_executable(TestErick
# #     meshVisualization.cpp
# #     DarcySpaces.cpp
# #     mainDarcy.cpp
# # )
# add_executable(TestErick mainDarcy.cpp meshVisualization.cpp)
# add_executable(TestDarcy DarcySpaces.cpp)

# target_link_libraries(TestDarcy NeoPZ::pz)
# target_link_libraries(TestErick NeoPZ::pz)

add_executable(TestErick meshVisualization.cpp) 
add_executable(TestDarcy DarcySpaces.cpp)
target_link_libraries(TestErick NeoPZ::pz)
target_link_libraries(TestDarcy NeoPZ::pz)



