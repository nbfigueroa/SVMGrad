cmake_minimum_required(VERSION 2.8.3)
project(svm_grad)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

find_package(catkin REQUIRED)
find_package(Armadillo REQUIRED)

## System dependencies are found with CMake's conventions
# find_package(Boost REQUIRED COMPONENTS system)

## Uncomment this if the package has a setup.py. This macro ensures
## modules and global scripts declared therein get installed
## See http://ros.org/doc/api/catkin/html/user_guide/setup_dot_py.html
# catkin_python_setup()

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need
catkin_package(
  INCLUDE_DIRS include/svm_grad
  LIBRARIES svm_grad
#  CATKIN_DEPENDS other_catkin_pkg
#  DEPENDS system_lib
)

###########
## Build ##
###########

## Specify additional locations of header files
## Your package locations should be listed before other locations
include_directories(
   include/svm_grad 
   ${ARMADILLO_INCLUDE_DIRS}
)

## Declare a C++ library
add_library(svm_grad src/svm_grad/svm_grad.cpp)
target_link_libraries (svm_grad 
	${catkin_LIBRARIES}
	armadillo)


## Declare a C++ executable
add_executable(test_svmgrad src/test_svmgrad.cpp)

## Specify libraries to link a library or executable target against
target_link_libraries(test_svmgrad
   svm_grad
   ${catkin_LIBRARIES}
 )


