cmake_minimum_required(VERSION 3.1)
set( CMAKE_CXX_COMPILER "g++" )
set( CMAKE_C_COMPILER "gcc" )
project(fcpw-libigl-fpdc-example)

list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#Enable optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


include(libigl)
add_subdirectory(fcpw)
add_executable(${PROJECT_NAME} main.cpp fpdc/polygon_wall.cpp)

target_link_libraries(${PROJECT_NAME} fcpw igl::core)
target_include_directories(${PROJECT_NAME} PUBLIC ${FCPW_EIGEN_INCLUDES})
target_include_directories(${PROJECT_NAME} PUBLIC ${FCPW_ENOKI_INCLUDES})
target_include_directories(${PROJECT_NAME} PUBLIC fpdc)


