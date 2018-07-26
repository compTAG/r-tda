cmake_minimum_required(VERSION 2.8)

find_package(Boost)

set(DIONYSUS_INCLUDE_DIR "../dionysus/include")
include_directories(${DIONYSUS_INCLUDE_DIR} ${Boost_INCLUDE_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(filtration-persistence  filtration-persistence.cpp)
target_link_libraries(filtration-persistence ${libraries})
