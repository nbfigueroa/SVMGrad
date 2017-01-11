cmake_minimum_required(VERSION 2.4.6)
PROJECT(ASVMLearning)
SET(PROJECT_NAME ASVMLearning)

set(USE_NLOPT CACHE BOOL ON)

if(USE_NLOPT)

	find_library(NLOPT_LIB NAMES nlopt nlopt_cxx libnlopt-0 PATHS "$ENV{PATH}" "$ENV{LD_LIBRARY_PATH}"  "$ENV{DYLD_LIBRARY_PATH}" DOC "NLOPT library") 
	
	if(NLOPT_LIB)
		MESSAGE(STATUS "NLOPT lib found ${NLOPT_LIB}")
		GET_FILENAME_COMPONENT( NLOPT_LIB_PATH ${NLOPT_LIB} PATH)
		GET_FILENAME_COMPONENT( NLOPT_LIB_NAME ${NLOPT_LIB} NAME_WE)
		STRING(REPLACE "lib" "" NLOPT_LIB_NAME ${NLOPT_LIB_NAME})
	else(NLOPT_LIB)
		MESSAGE(ERROR "\nCannot find libnlopt.so in standard locations. Please\n - add the containing folder to your LD_LIBRARY_PATH (linux); PATH (windows);DYLD_LIBRARY_PATH(mac)   \n OR \n - set USE_NLOPT to False")
	endif(NLOPT_LIB)
	
	FIND_PATH(NLOPT_INCLUDE_PATH nlopt.hpp PATHS "$ENV{PATH}" DOC "Folder containing the header file nlopt.hpp")
	if(NLOPT_INCLUDE_PATH)
		MESSAGE(STATUS "nlopt.hpp found in ${NLOPT_INCLUDE_PATH}")
	else(NLOPT_INCLUDE_PATH)
		MESSAGE(ERROR "\nCannot find nlopt.hpp in standard locations. Please\n - Set NLOPT_INCLUDE_PATH manually \n OR \n - Set USE_NLOPT to OFF")
	endif(NLOPT_INCLUDE_PATH)
endif(USE_NLOPT)



#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin )
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib )

if(UNIX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-result")
endif(UNIX)

set(SRC_FILES src/asvm.cpp  src/asvmdata.cpp  src/asvm_smo_solver.cpp src/util.cpp src/svm.cpp)
if(USE_NLOPT)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_NLOPT")
	list(APPEND SRC_FILES src/asvm_nlopt_solver.cpp src/quadprog_nlopt_solver.cpp)
endif(USE_NLOPT)

#common commands for building c++ executables and libraries
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include )
if(USE_NLOPT)
	INCLUDE_DIRECTORIES(${NLOPT_INCLUDE_PATH})
endif(USE_NLOPT)

if(WIN32)
add_library(${PROJECT_NAME}   ${SRC_FILES})
else(WIN32)
add_library(${PROJECT_NAME} SHARED ${SRC_FILES})
endif(WIN32)


if(USE_NLOPT)
  target_link_libraries(${PROJECT_NAME} ${NLOPT_LIB} )
endif(USE_NLOPT)

add_executable(train src/main.cpp)
if(USE_NLOPT)
  target_link_libraries(train  ${PROJECT_NAME} ${NLOPT_LIB} )
else(USE_NLOPT)
	target_link_libraries(train ${PROJECT_NAME} )
endif(USE_NLOPT)


# Generating Config.cmake file for "Find"-ing from other cmake projects. 
SET(INC_DIRS ${NLOPT_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/include)
SET(LIB_DIRS ${PROJECT_SOURCE_DIR}/lib ${NLOPT_LIB_PATH})
SET(LIBS ${NLOPT_LIB_NAME} ${PROJECT_NAME})
if(USE_NLOPT)
SET(NLOPT_CHOICE true)
else(USE_NLOPT)
SET(NLOPT_CHOICE false)
endif(USE_NLOPT)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}.in ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}Config.cmake @ONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/matlab/SVMWidget/make.in ${CMAKE_CURRENT_SOURCE_DIR}/matlab/SVMWidget/make.m @ONLY)

if(UNIX)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/run_matlab.in ${CMAKE_CURRENT_SOURCE_DIR}/run_matlab.sh @ONLY)
EXECUTE_PROCESS(COMMAND chmod a+x run_matlab.sh WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif(UNIX)