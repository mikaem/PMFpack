cmake_minimum_required(VERSION 2.6)

find_package(Boost COMPONENTS system filesystem timer REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})
#------------------------------------------------------------------------------

list(APPEND CXX_DEFINITIONS "-O1 -std=c++11 -Wno-narrowing")

add_definitions(${CXX_DEFINITIONS})
add_executable(demo test_multiprecision.cpp)
target_link_libraries(demo ${Boost_LIBRARIES} )
