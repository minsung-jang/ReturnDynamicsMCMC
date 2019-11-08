cmake_minimum_required(VERSION 2.8.3)
project(ReturnDynamicsMCMC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") # Add -std=c++11

## Eigen Package
find_package(Eigen3 REQUIRED)

## Include directories for header files
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB CPPS "src/*.cpp")
file(GLOB HEADERS "/include/*.h")

#Otherwise, directly pick up certain names of source files
set(SOURCE_FILES
    ${CPPS}
    ${HEADERS})

## Add Library

## Target Library

## Executable files
add_executable(${PROJECT_NAME} ${SOURCE_FILES})
