# 
# study the pnp method using 2d and 3d points  
#

project(pnp_hybrid)
cmake_minimum_required(VERSION 2.8.3)

set(CMAKE_CXX_FLAGS "-ggdb -g -O3 -w -std=c++11 -fpermissive") # -fPIC -fpermissive

find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS} "/usr/include/eigen3" "./src")
add_subdirectory(src bin)
