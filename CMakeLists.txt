# 
# study the pnp method using 2d and 3d points  
#

project(pnp_hybrid)
cmake_minimum_required(VERSION 2.8.3)

set(CMAKE_CXX_FLAGS "-ggdb -g -O3 -w -std=c++11 -fpermissive") # -fPIC -fpermissive

##### ceres #####
find_package(Ceres REQUIRED)
include_directories(${CERES_INCLUDE_DIRS} "/usr/include/eigen3" "./src")

##### opencv ##### 
find_package(OpenCV REQUIRED)

##### sophus #####
find_package( sophus REQUIRED )
include_directories( ${sophus_INCLUDE_DIRS} )

add_subdirectory(src bin)
