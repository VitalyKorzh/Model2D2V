cmake_minimum_required(VERSION 3.27.4 FATAL_ERROR)

project(model)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

find_package(LinSolver REQUIRED)
#find_package(OpenMP REQUIRED)

file(GLOB SRC ${CMAKE_SOURCE_DIR}/src/*.cpp)

add_executable(${PROJECT_NAME} ${SRC})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include
${CMAKE_SOURCE_DIR}/include/linear_math)
#target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(${PROJECT_NAME} PRIVATE LinSolver::lin_solver 
                                                )