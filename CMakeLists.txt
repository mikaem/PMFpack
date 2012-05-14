cmake_minimum_required(VERSION 2.6)

project(PMF)
set(PMF_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include" CACHE PATH "C/C++ header installation directory.")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR} CACHE PATH "PMF install prefix" FORCE
)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

#------------------------------------------------------------------------------
# Install header files

file(GLOB _HEADERS ${PMF_INCLUDE_DIR}/*.h)
list(APPEND HEADERS ${_HEADERS})
install(FILES ${_HEADERS} DESTINATION include COMPONENT Development)
include_directories(${PMF_INCLUDE_DIR})

#------------------------------------------------------------------------------
# Find sources
file(GLOB CPP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
list(APPEND SOURCES ${CPP_SOURCES})
file(GLOB C_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)
list(APPEND SOURCES ${C_SOURCES})

add_definitions("-fPIC")
# Dependency
find_package(GSL)
include_directories(${GSL_INCLUDE_DIRS})
add_library(PMF ${HEADERS} ${SOURCES})
string(STRIP GSL_LIBRARIES "${GSL_LIBRARIES}")
target_link_libraries(PMF ${GSL_LIBRARIES})
install(TARGETS PMF DESTINATION lib)

add_subdirectory(demo)

add_subdirectory(swig)
