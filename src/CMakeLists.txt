cmake_minimum_required(VERSION 3.12)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(Boost_INCLUDE_DIR /usr/include)

find_package(Boost 1.53)
project(dragom)

if(Boost_FOUND)
  add_executable(
    dragom.exe
    main.cpp
    SG.cpp
    getSG.cpp
    mergeSG.cpp
    extendAnchor.cpp
    unionFind.cpp
  )
else()
  message(FATAL_ERROR "You must have boost library (>=1.53) to compile DRAGoM!")
endif()
