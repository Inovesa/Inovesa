project(inovesa)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(./src SRC_LIST)

set(SRC_LIST
			./src/CL/CLProgApplyHM.cpp
			./src/CL/CLProgs.cpp
			./src/CL/OpenCLHandler.cpp
			./src/HM/HeritageMap.cpp
			./src/HM/FokkerPlanckMap.cpp
			./src/HM/RotationMap.cpp
			./src/HM/WakeKickMap.cpp
			./src/IO/HDF5File.cpp
			./src/Display.cpp
			./src/PhaseSpace.cpp
			./src/main.cpp
			./src/Share.cpp
)

set(INC_LIST
			./inc/CL/CLProgApplyHM.hpp
			./inc/CL/CLProgs.hpp
			./inc/CL/OpenCLHandler.hpp
			./inc/HM/HeritageMap.hpp
			./inc/HM/FokkerPlanckMap.hpp
			./inc/HM/Identity.hpp
			./inc/HM/RotationMap.hpp
			./inc/HM/WakeKickMap.hpp
			./inc/IO/HDF5File.hpp
			./inc/Display.hpp
			./inc/PhaseSpace.hpp
			./inc/Ruler.hpp
			./inc/Share.hpp
			./inc/defines.hpp
			./inc/fixed_point.h
)
add_executable(${PROJECT_NAME} ${SRC_LIST} ${INC_LIST})

add_definitions(-std=c++11)
if(DEFINED CMAKE_BUILD_TYPE)
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
		add_definitions(-O3 -Wall -Wextra )
        MESSAGE ("Doing release compile, using optimization flags.")
    else()
		add_definitions(-g -Wall -Wextra -DDEBUG)
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
find_package(HDF5 COMPONENTS CXX REQUIRED)
find_package(PNG REQUIRED)
include_directories(./inc
	${GLEW_INCLUDE_DIRS} ${GLFW_INCLUDE_DIRS} ${GLM_INCLUDE_DIRS}
	${HDF5_INCLUDE_DIRS}
	${OPENCL_INCLUDE_DIRS}
)
SET(LIBS ${LIBS} ${GSL_LIBRARIES} ${Boost_LIBRARIES}
	${GLFW_LIBRARIES} ${GLEW_LIBRARIES} ${GL_LIBRARY}
	${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES}
	${OPENGL_LIBRARIES} ${OPENCL_LIBRARIES}
	${PNG_LIBRARIES}
)
target_link_libraries(${PROJECT_NAME} ${LIBS})
