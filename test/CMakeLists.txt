cmake_minimum_required(VERSION 2.6)

add_executable(demo demo.cpp)

include_directories("${PROJECT_BINARY_DIR}/include")
include_directories(${GSL_INCLUDE_DIRS})
target_link_libraries(demo_root ${GSL_LIBRARIES})
target_link_libraries(demo_root PMF)
