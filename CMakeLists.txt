cmake_minimum_required(VERSION 3.1.0)
project (fdtd1d)

#set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

include_directories(include)

set (PROJECT_LINK_LIBS "${PROJECT_LINK_LIBS} " )
link_directories(/usr/lib)

file(GLOB SOURCES "src/*.cc")
add_executable(fdtd1d ${SOURCES})

target_link_libraries(fdtd1d -lpthread -fopenmp)


