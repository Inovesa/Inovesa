project(inovesa)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(./src SRC_LIST)

set(SRC_LIST
            ./src/CL/OpenCLHandler.cpp
            ./src/HM/HeritageMap.cpp
            ./src/HM/FokkerPlanckMap.cpp
            ./src/HM/KickMap.cpp
            ./src/HM/RFKickMap.cpp
            ./src/HM/RotationMap.cpp
            ./src/HM/WakeKickMap.cpp
            ./src/HM/WakePotentialMap.cpp
            ./src/HM/WakeFunctionMap.cpp
            ./src/IO/Display.cpp
            ./src/IO/GUI/GUIElement.cpp
            ./src/IO/GUI/Plot2DLine.cpp
            ./src/IO/GUI/Plot3DColormap.cpp
            ./src/IO/HDF5File.cpp
            ./src/IO/ProgramOptions.cpp
            ./src/ElectricField.cpp
            ./src/Impedance.cpp
            ./src/PhaseSpace.cpp
            ./src/main.cpp
)

set(INC_LIST
            ./inc/CL/local_cl.hpp
            ./inc/CL/OpenCLHandler.hpp
            ./inc/HM/HeritageMap.hpp
            ./inc/HM/FokkerPlanckMap.hpp
            ./inc/HM/KickMap.hpp
            ./inc/HM/RFKickMap.hpp
            ./inc/HM/Identity.hpp
            ./inc/HM/RotationMap.hpp
            ./inc/HM/WakePotentialMap.hpp
            ./inc/HM/WakeKickMap.hpp
            ./inc/HM/WakeFunctionMap.hpp
            ./inc/IO/Display.hpp
            ./inc/IO/GUI/GUIElement.hpp
            ./inc/IO/GUI/Plot2DLine.hpp
            ./inc/IO/GUI/Plot3DColormap.hpp
            ./inc/IO/HDF5File.hpp
            ./inc/IO/ProgramOptions.hpp
            ./inc/ElectricField.hpp
            ./inc/Impedance.hpp
            ./inc/PhaseSpace.hpp
            ./inc/Ruler.hpp
            ./inc/defines.hpp
            ./inc/fixed_point.h
)
set(MISC_LIST
    default.cfg
)
add_executable(${PROJECT_NAME} ${SRC_LIST} ${INC_LIST} ${MISC_LIST})


execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions("-DGIT_BRANCH=\"${GIT_BRANCH}\"")

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
find_package(OpenCL)
IF(OPENCL_FOUND)
    add_definitions( -DINOVESA_USE_CL)
    MESSAGE ("Found OpenCL. Will add support.")
    find_package(clFFT)
    IF(clFFT_FOUND)
        add_definitions( -DINOVESA_USE_CLFFT)
        MESSAGE ("Found clFFT. Will add support.")
    ENDIF()
ENDIF()
find_package(OpenGL)
IF(OPENGL_FOUND)
    add_definitions( -DINOVESA_USE_GUI)
    MESSAGE ("Found OpenGL. Will add support.")
    find_package(GLEW REQUIRED)
    find_package(GLFW REQUIRED)
    add_definitions( -DGLFW_VERSION_MAJOR=${GLFW_VERSION_MAJOR})
ENDIF()
find_package(Boost COMPONENTS filesystem program_options system REQUIRED)
find_package(FFTW REQUIRED)
find_package(HDF5 COMPONENTS CXX)
IF(HDF5_FOUND)
    add_definitions( -DINOVESA_USE_HDF5)
    MESSAGE ("Found HDF5. Will add support.")
    set (INOVESA_HDF5_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
    set (INOVESA_HDF5_LIBRARIES ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
ENDIF()
find_package(PNG)
IF(PNG_FOUND)
    add_definitions( -DINOVESA_USE_PNG)
    MESSAGE ("Found PNG. Will add support.")
ENDIF()

include_directories(./inc
    ${GLEW_INCLUDE_DIRS} ${GLFW_INCLUDE_DIRS} ${GLM_INCLUDE_DIRS}
    ${INOVESA_HDF5_INCLUDE_DIRS} ${PNG_INCLUDE_DIRS} ${FFTW_INCLUDE_DIRS}
    ${OPENCL_INCLUDE_DIRS} ${CLFFT_INCLUDE_DIRS}
)
SET(LIBS ${LIBS} ${Boost_LIBRARIES}
    ${GLFW_LIBRARIES} ${GLEW_LIBRARIES}
    ${GL_LIBRARY} ${INOVESA_HDF5_LIBRARIES}
    ${OPENGL_LIBRARIES} ${OPENCL_LIBRARIES}
    ${PNG_LIBRARIES} ${FFTW_LIBRARIES}
    ${CLFFT_LIBRARIES}
)
target_link_libraries(${PROJECT_NAME} ${LIBS})
