cmake_minimum_required(VERSION 2.6)

enable_language (Fortran)

include_directories(${PMF_INCLUDE_DIR})

add_executable(driver DRIVER.f)

target_link_libraries(driver fortPMF)

