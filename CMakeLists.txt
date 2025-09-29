cmake_minimum_required(VERSION 3.16)

project(FluidTut LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(SOURCES 
      src/main.cpp
      src/cell.cpp
      src/grid2D.cpp
      src/isentropicVortex.cpp
      src/solver2D.cpp
      src/rungeKutta3.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})

target_include_directories(${PROJECT_NAME}
   PRIVATE
      ${CMAKE_SOURCE_DIR}/include)


