cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(HololensIGTLink)

#SET(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/../Slicer-unixBuild/VTKv6-build")
find_package(VTK 7.0 REQUIRED)
include(${VTK_USE_FILE})

SET(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/../../OpenIGTLinkPointCloud/OpenIGTLinkPointCloud-Build")
message(${CMAKE_PREFIX_PATH})
find_package(OpenIGTLink REQUIRED)
include(${OpenIGTLink_USE_FILE})

add_executable (PolyDataReceiver PolyDataReceiver.cpp)
target_link_libraries (PolyDataReceiver OpenIGTLink ${VTK_LIBRARIES})

