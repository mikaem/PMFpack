cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/modules")
find_package(Boost COMPONENTS system filesystem timer REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# find_package(FFTW)
# if (FFTW_FOUND)
#     include_directories(${FFTW_INCLUDE_DIR})
# endif()

#------------------------------------------------------------------------------

list(APPEND CXX_DEFINITIONS "-O1 -std=c++11 -Wno-narrowing")

add_definitions(${CXX_DEFINITIONS})
add_executable(demo test_multiprecision.cpp)
target_link_libraries(demo ${Boost_LIBRARIES} )
# if (FFTW_FOUND)
#     target_link_libraries(demo ${FFTW_LIBRARY} ${FFTW_THREADS_LIBRARY})
# endif()
