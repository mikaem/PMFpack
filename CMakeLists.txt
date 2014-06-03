cmake_minimum_required(VERSION 2.6)

project(PMF)
set(PMF_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include" CACHE PATH "C/C++ header installation directory.")
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR} CACHE PATH "PMF install prefix" FORCE
)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# Enable Boost?
set(PMF_ENABLE_BOOST ${PMF_ENABLE_BOOST} CACHE BOOL "FALSE")

find_package(Boost COMPONENTS system filesystem timer REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

#------------------------------------------------------------------------------
# Install header files
file(GLOB _HEADERS ${PMF_INCLUDE_DIR}/*.h)
list(APPEND HEADERS ${_HEADERS})
if (NOT PMF_ENABLE_BOOST)
  file(GLOB Boost_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/Boost*.h)
  list(REMOVE_ITEM HEADERS ${Boost_HEADERS})  
endif()
install(FILES ${HEADERS} DESTINATION include COMPONENT Development)
include_directories(${PMF_INCLUDE_DIR})
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/multiprecision")
#------------------------------------------------------------------------------
# Find sources
file(GLOB CPP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
list(APPEND SOURCES ${CPP_SOURCES})
file(GLOB C_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)
list(APPEND SOURCES ${C_SOURCES})
if (NOT PMF_ENABLE_BOOST)
  file(GLOB Boost_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/Boost*.cpp)
  list(REMOVE_ITEM SOURCES ${Boost_SOURCES})  
endif()

list(APPEND PMF_CXX_DEFINITIONS "-fPIC -O1 -std=c++11 -Wno-narrowing")

if (PMF_ENABLE_BOOST)
  list(APPEND PMF_CXX_DEFINITIONS "-DHAS_BOOST")
  if (CMAKE_BUILD_TYPE STREQUAL "Developer" OR CMAKE_BUILD_TYPE STREQUAL "Debug")
    # Makes it print out a lot of debug info on Boost:
    list(APPEND PMF_CXX_DEFINITIONS "-DBOOST_MATH_INSTRUMENT")
  endif()
endif()


# Dependency
Find_package(GSL)
add_definitions(${PMF_CXX_DEFINITIONS})
include_directories(${GSL_INCLUDE_DIRS})
add_library(PMF ${HEADERS} ${SOURCES})
string(STRIP GSL_LIBRARIES "${GSL_LIBRARIES}")
target_link_libraries(PMF ${GSL_LIBRARIES})

#FIND_PACKAGE (HDF5 NAMES hdf5)
FIND_PACKAGE (HDF5 HINTS ${HDF5_DIR}) # Find non-cmake built HDF5
INCLUDE_DIRECTORIES (${HDF5_INCLUDE_DIR})
target_link_libraries( PMF ${LINK_LIBS} ${HDF5_LIBRARIES} z ${Boost_LIBRARIES})

install(TARGETS PMF DESTINATION lib)

add_subdirectory(demo)

add_subdirectory(fortran)

add_subdirectory(demofortran)

add_subdirectory(swig)

