cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(Compare)

set (CMAKE_PREIX_PATH /home/labuser/libtorch)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -fpermissive")
set(CMAKE_CXX_FLAGS_DEBUG "-g -pg")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_RL_FLAGS} ")

add_executable(generate src/Generate.cpp
                    include/Wavepacket.hpp
                    include/Environment_MhollandRL.hpp
                    )

target_include_directories(generate PRIVATE include)

target_link_libraries(generate  Eigen3::Eigen)

set_property(TARGET generate  PROPERTY CXX_STANDARD 23)
