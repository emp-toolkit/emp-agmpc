cmake_minimum_required (VERSION 2.8.11)
project (emp-agmpc)
set(NAME "emp-agmpc")

find_path(CMAKE_FOLDER NAMES cmake/emp-tool-config.cmake)
include(${CMAKE_FOLDER}/cmake/common.cmake)
include(${CMAKE_FOLDER}/cmake/enable_rdseed.cmake)
include(${CMAKE_FOLDER}/cmake/enable_float.cmake)

FIND_PACKAGE(emp-ot REQUIRED)
INCLUDE_DIRECTORIES(${EMP-OT_INCLUDE_DIRS})

install(FILES cmake/emp-agmpc-config.cmake DESTINATION cmake/)
install(DIRECTORY emp-agmpc DESTINATION include/)

ENABLE_TESTING()
ADD_SUBDIRECTORY(test)
