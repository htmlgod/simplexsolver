cmake_minimum_required(VERSION 3.4)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(simplexsolver)

add_library(tools STATIC
        ${CMAKE_CURRENT_SOURCE_DIR}/sources/tools.cpp
)

add_library(simplex STATIC
        ${CMAKE_CURRENT_SOURCE_DIR}/sources/simplex.cpp
)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(simplexsolver
        ${CMAKE_CURRENT_SOURCE_DIR}/sources/main.cpp
)

target_link_libraries(simplexsolver simplex tools)
