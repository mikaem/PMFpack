cmake_minimum_required(VERSION 2.6)

add_executable(demo demo.cpp)

include_directories("${PROJECT_BINARY_DIR}/include" ${GSL_INCLUDE_DIRS})
target_link_libraries(demo ${GSL_LIBRARIES} PMF ${Boost_LIBRARIES} )
