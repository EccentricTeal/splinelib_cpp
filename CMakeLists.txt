cmake_minimum_required(VERSION 3.0.2)
project(splinelib_cpp)
add_compile_options(-std=c++14)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
)

catkin_package(
  INCLUDE_DIRS include
  LIBRARIES nurbslib
#  CATKIN_DEPENDS roscpp rospy std_msgs
#  DEPENDS eigen
)

include_directories(
  include
  ${catkin_INCLUDE_DIRS}
  /usr/include/eigen3/
)

## Declare a C++ library
# add_library(${PROJECT_NAME}
#   src/${PROJECT_NAME}/routegen_test.cpp
# )

add_library(nurbslib
  src/nurbs.cc
)
target_link_libraries(nurbslib
  ${catkin_LIBRARIES}
)

