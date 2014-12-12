project(fast_rotate)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(./src SRC_LIST)

set(SRC_LIST
			./src/CL/CLProgApplyHM.cpp
			./src/CL/CLProgRotateKick.cpp
			./src/CL/CLProgs.cpp
			./src/CL/OpenCLHandler.cpp
			./src/Display.cpp
			./src/main.cpp
)

set(INC_LIST
			./inc/CL/CLProgApplyHM.hpp
			./inc/CL/CLProgRotateKick.hpp
			./inc/CL/CLProgs.hpp
			./inc/CL/OpenCLHandler.hpp
			./inc/Display.hpp
			./inc/Mesh2D.hpp
			./inc/Ruler.hpp
)
add_executable(${PROJECT_NAME} ${SRC_LIST} ${INC_LIST})

add_definitions(-std=c++11)
if(DEFINED CMAKE_BUILD_TYPE)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        add_definitions(-O3)
        MESSAGE ("Doing release compile, using optimization flags.")
    else()
        add_definitions(-g -Wall -DDEBUG)
        MESSAGE ("Doing debug compile, including debug symbols.")
    endif()
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "./cmake")
find_package(OpenCL REQUIRED)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(GLFW REQUIRED)
find_package(GSL REQUIRED)
find_package(GSL COMPONENTS linalg REQUIRED)
find_package(Boost COMPONENTS filesystem system REQUIRED)
include_directories(./inc ${GLEW_INCLUDE_DIRS} ${GLFW_INCLUDE_DIRS} ${GLM_INCLUDE_DIRS} ${OPENCL_INCLUDE_DIRS})
SET(LIBS ${LIBS} ${GSL_LIBRARIES} ${Boost_LIBRARIES} ${GLFW_LIBRARIES} ${GLEW_LIBRARIES} ${GL_LIBRARY} ${OPENGL_LIBRARIES} ${OPENCL_LIBRARIES})
target_link_libraries(${PROJECT_NAME} ${LIBS})
