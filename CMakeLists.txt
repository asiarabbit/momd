######################################################
# SUNNY project, Anyang Normal University, IMP-CAS
#
# CMakeLists.txt for user-composed programs
# Author: SUN Yazhou, asia.rabbit@163.com
# Created: 2020/07/08
# Last modified: 2020/07/08, SUN Yazhou
######################################################

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
set(PROJECT_NAME MOMD)
project(${PROJECT_NAME})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules")
set(CMAKE_BUILD_TYPE "Debug")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -std=c++11 -O0 -Wall -g -ggdb")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

add_subdirectory(momd) # library path
add_subdirectory(src)   # user-defined source file path
