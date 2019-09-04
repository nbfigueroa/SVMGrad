cmake_minimum_required(VERSION 3.9)

# Options. Turn on with 'cmake -Druntests=ON'.
option(runtests "Build all tests." OFF) # Makes boolean 'test' available.

project(svm_grad)

# Default to C99
if(NOT CMAKE_C_STANDARD)
  	set(CMAKE_C_STANDARD 99)
endif()

# Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  	set(CMAKE_CXX_STANDARD 14)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  	add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(Armadillo REQUIRED)

include_directories(
   include
   ${ARMADILLO_INCLUDE_DIRS}
)

set(CORE_SOURCES
	src/svm_grad.cpp
)

add_library(${PROJECT_NAME} SHARED 
  ${CORE_SOURCES}
)

target_link_libraries(${PROJECT_NAME}
	armadillo
)

install(DIRECTORY include/
  	DESTINATION include
)

install(TARGETS ${PROJECT_NAME}
  	ARCHIVE DESTINATION lib
  	LIBRARY DESTINATION lib
  	RUNTIME DESTINATION bin
)

if (runtests)
    if (APPLE)
      add_definitions(-DGTEST_USE_OWN_TR1_TUPLE)
      add_definitions(-D__GLIBCXX__)
    endif (APPLE)

    enable_testing()

    add_executable(test_svmgrad tests/test_svmgrad.cpp)
    target_link_libraries(test_svmgrad
      ${PROJECT_NAME}
    )
    add_test(NAME test_svmgrad COMMAND test_svmgrad)
endif()

