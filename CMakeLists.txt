cmake_minimum_required(VERSION 3.10)

project(inovesa)

## User configuration
# Enable interrupts (Ctrl+C)
add_definitions(-DINOVESA_ENABLE_INTERRUPT=1)


# (default) options for Debug and Release builds
SET( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${-Wall}" )
if((DEFINED CMAKE_BUILD_TYPE) AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
  if (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
    add_definitions(-Wextra)
	add_definitions(-g)
  endif()
  add_definitions(-DDEBUG=1)
  MESSAGE ("Doing debug compile, including debug symbols.")
  unset(QUIET)
else()
  if (("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
    add_definitions(-O3 )
  endif()
  add_definitions(-w )
  MESSAGE ("Doing release compile, using optimization flags.")
  set(QUIET QUIET)
endif()

## end of user configuration

set (INOVESA_VERSION_MAJOR 1)
set (INOVESA_VERSION_MINOR 2)
set (INOVESA_VERSION_FIX -1)
# pre release numbers:
# -1: alpha
# -2: beta
# -3 RC1
# -n RC(n-2)

SET (COVERAGE OFF CACHE BOOL "Coverage")

configure_file (
  "${PROJECT_SOURCE_DIR}/InovesaConfig.hpp.in"
  "${PROJECT_BINARY_DIR}/InovesaConfig.hpp"
)
include_directories("${PROJECT_BINARY_DIR}")

aux_source_directory(./src SRC_LIST)

set(SRC_LIST
  ./src/CL/CLProfiler.cpp
  ./src/CL/OpenCLHandler.cpp
  ./src/IO/Display.cpp
  ./src/IO/FSPath.cpp
  ./src/IO/GUI/GUIElement.cpp
  ./src/IO/GUI/Plot1DLine.cpp
  ./src/IO/GUI/Plot2DPoints.cpp
  ./src/IO/GUI/Plot2DPrimitive.cpp
  ./src/IO/GUI/Plot3DColormap.cpp
  ./src/IO/HDF5File.cpp
  ./src/IO/ProgramOptions.cpp
  ./src/PS/ElectricField.cpp
  ./src/PS/PhaseSpace.cpp
  ./src/PS/PhaseSpaceFactory.cpp
  ./src/SM/DriftMap.cpp
  ./src/SM/FokkerPlanckMap.cpp
  ./src/SM/KickMap.cpp
  ./src/SM/DynamicRFKickMap.cpp
  ./src/SM/RFKickMap.cpp
  ./src/SM/RotationMap.cpp
  ./src/SM/Identity.cpp
  ./src/SM/SourceMap.cpp
  ./src/SM/WakeKickMap.cpp
  ./src/SM/WakePotentialMap.cpp
  ./src/Z/CollimatorImpedance.cpp
  ./src/Z/ConstImpedance.cpp
  ./src/Z/Impedance.cpp
  ./src/Z/ImpedanceFactory.cpp
  ./src/Z/FreeSpaceCSR.cpp
  ./src/Z/ParallelPlatesCSR.cpp
  ./src/Z/ResistiveWall.cpp
  ./src/FFTWWrapper.cpp
  ./src/HelperFunctions.cpp
)

set(INC_LIST
  ./inc/CL/CLProfiler.hpp
  ./inc/CL/local_cl.hpp
  ./inc/CL/OpenCLHandler.hpp
  ./inc/PS/ElectricField.hpp
  ./inc/PS/PhaseSpace.hpp
  ./inc/PS/PhaseSpaceFactory.hpp
  ./inc/PS/Ruler.hpp
  ./inc/IO/Display.hpp
  ./inc/IO/FSPath.hpp
  ./inc/IO/GUI/GUIElement.hpp
  ./inc/IO/GUI/Plot1DLine.hpp
  ./inc/IO/GUI/Plot2DPoints.hpp
  ./inc/IO/GUI/Plot2DPrimitive.hpp
  ./inc/IO/GUI/Plot3DColormap.hpp
  ./inc/IO/HDF5File.hpp
  ./inc/IO/ProgramOptions.hpp
  ./inc/SM/DriftMap.hpp
  ./inc/SM/FokkerPlanckMap.hpp
  ./inc/SM/KickMap.hpp
  ./inc/SM/DynamicRFKickMap.hpp
  ./inc/SM/RFKickMap.hpp
  ./inc/SM/RotationMap.hpp
  ./inc/SM/Identity.hpp
  ./inc/SM/SourceMap.hpp
  ./inc/SM/WakePotentialMap.hpp
  ./inc/SM/WakeKickMap.hpp
  ./inc/Z/CollimatorImpedance.hpp
  ./inc/Z/ConstImpedance.hpp
  ./inc/Z/Impedance.hpp
  ./inc/Z/ImpedanceFactory.hpp
  ./inc/Z/FreeSpaceCSR.hpp
  ./inc/Z/ParallelPlatesCSR.hpp
  ./inc/Z/ResistiveWall.hpp
  ./inc/FFTWWrapper.hpp
  ./inc/HelperFunctions.hpp
  ./inc/defines.hpp
)

add_executable(${PROJECT_NAME} ./src/main.cpp ${SRC_LIST} ${INC_LIST} ${MISC_LIST})

include(CheckIncludeFileCXX)

execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_BRANCH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions("-DGIT_BRANCH=\"${GIT_BRANCH}\"")

execute_process(
  COMMAND git rev-parse --verify HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
add_definitions("-DGIT_COMMIT=\"${GIT_COMMIT}\"")
include_directories(./inc)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# workaround CMake bug where previous definition is ignored
set(CMAKE_CXX_FLAGS "-std=c++14 ${CMAKE_CXX_FLAGS}")

# make compatible to newer boost and cmake
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 NEW)
endif()

# workaround problem in boost with g++
IF( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
    set(CMAKE_CXX_FLAGS "-fext-numeric-literals ${CMAKE_CXX_FLAGS}")
ENDIF()

# turn on automatic optimizations for current architecture when GCC is used
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  add_definitions(-march=native)
  MESSAGE ("Builing for current architecture. Binary might not work on other CPU.")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


## OpenCL (optional)
IF((NOT DEFINED ENV{INOVESA_USE_OPENCL}) OR ("$ENV{INOVESA_USE_OPENCL}" EQUAL 1))
find_package(OpenCL ${QUIET})
IF(OPENCL_FOUND)
    add_definitions( -DINOVESA_USE_OPENCL=1)
    MESSAGE ("Found OpenCL. Will add support.")
    SET(LIBS ${LIBS} ${OPENCL_LIBRARIES})
    include_directories(${OPENCL_INCLUDE_DIRS})
    find_package(clFFT ${QUIET})
    IF(clFFT_FOUND)
        add_definitions( -DINOVESA_USE_CLFFT=1)
        MESSAGE ("Found clFFT. Will add support.")
        SET(LIBS ${LIBS} ${CLFFT_LIBRARIES})
        include_directories(${CLFFT_INCLUDE_DIRS})
    ELSE()
        add_definitions( -DINOVESA_USE_CLFFT=0)
    ENDIF()
    IF((NOT DEFINED INOVESA_ENABLE_CLPROFILING))
        add_definitions(-DINOVESA_ENABLE_CLPROFILING=0)
    ENDIF()
ELSE()
    add_definitions( -DINOVESA_USE_OPENCL=0)
    MESSAGE ("Did not find OpenCL. Will compile without OpenCL support.")
ENDIF()
ELSE()
    MESSAGE ("Excluded OpenCL. Will compile without OpenCL support.")
ENDIF()

## OpenGL (optional)
IF((NOT DEFINED ENV{INOVESA_USE_OPENGL}) OR ("$ENV{INOVESA_USE_OPENGL}" EQUAL 1))
set(OpenGL_GL_PREFERENCE GLVND)
find_package(OpenGL ${QUIET})
find_package(GLEW ${QUIET})
find_package(GLFW ${QUIET})
IF((OPENGL_FOUND) AND (GLEW_FOUND) AND (GLFW_FOUND))
    add_definitions( -DINOVESA_USE_OPENGL=1)
    MESSAGE ("Found OpenGL. Will add support.")
    add_definitions( -DGLFW_VERSION_MAJOR=${GLFW_VERSION_MAJOR})
    include_directories(${GLEW_INCLUDE_DIRS} ${GLFW_INCLUDE_DIRS})
    set (LIBS ${LIBS} ${GLFW_LIBRARIES} ${GLEW_LIBRARIES}
              ${GL_LIBRARY} ${OPENGL_LIBRARIES}
    )
ELSE()
    add_definitions(-DINOVESA_USE_OPENGL=0)
    MESSAGE ("Did not find OpenGL. Will compile without OpenGL support.")
ENDIF()
ELSE()
    MESSAGE ("Excluded OpenGL. Will compile without OpenGL support.")
ENDIF()

## Boost (needed)
find_package(Boost COMPONENTS filesystem program_options system REQUIRED ${QUIET})
include_directories(${Boost_INCLUDE_DIRS})
SET(LIBS ${LIBS} ${Boost_LIBRARIES})

## FFTW (needed)
if (WIN32)
	find_package(FFTW3f CONFIG REQUIRED ${QUIET})
	SET(LIBS ${LIBS} FFTW3::fftw3f ${FFTW_LIBRARIES})
ELSE()
	find_package(FFTW REQUIRED ${QUIET})
	SET(LIBS ${LIBS} ${FFTW_LIBRARIES} )
ENDIF()
include_directories(${FFTW_INCLUDE_DIRS})

## HDF5 (optional)
find_package(HDF5 COMPONENTS CXX ${QUIET})
IF(HDF5_FOUND)
    add_definitions( -DINOVESA_USE_HDF5=1)
    MESSAGE ("Found HDF5. Will add support.")
    include_directories(${HDF5_INCLUDE_DIRS})
    set (LIBS ${LIBS} ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
ELSE()
    add_definitions( -DINOVESA_USE_HDF5=0)
    MESSAGE ("Did not find HDF5. Will compile without HDF5 support.")
ENDIF()

## PNG (optional)
find_package(OpenImageIO ${QUIET})
IF(OPENIMAGEIO_FOUND)
    add_definitions( -DINOVESA_USE_PNG=1)
    include_directories(${OPENIMAGEIO_INCLUDE_DIR})
    SET(LIBS ${LIBS} ${OPENIMAGEIO_LIBRARIES})
    MESSAGE ("Found OpenImageIO. Will add image I/O support.")
ELSE()
    add_definitions( -DINOVESA_USE_PNG=0)
    MESSAGE ("Did not find OpenImageIO. Will compile without PNG support.")
ENDIF()

target_link_libraries(${PROJECT_NAME} PRIVATE ${LIBS})

INSTALL(TARGETS inovesa RUNTIME DESTINATION bin)

find_package(Boost COMPONENTS unit_test_framework ${QUIET})
IF(Boost_UNIT_TEST_FRAMEWORK_FOUND)

    #Setup CMake to run tests
    enable_testing()
    include_directories (${Boost_INCLUDE_DIRS})

    #I like to keep test files in a separate source directory called test
    file(GLOB TEST_SRCS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} test/*.cpp)


    set(testName inovesa-test)

    #Add compile target
    add_executable(${testName} ${TEST_SRCS} ${SRC_LIST} ${INC_LIST} ${MISC_LIST})

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test/impedance.dat
        ${CMAKE_CURRENT_BINARY_DIR}/impedance.dat COPYONLY)

    #link to Boost libraries AND your targets and dependencies
    IF(COVERAGE AND ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU"))
        MESSAGE ("Build for coverage information.")
        target_compile_options(${testName} PRIVATE --coverage)
        target_link_libraries(${testName} PRIVATE ${Boost_LIBRARIES} ${LIBS} --coverage)
    ELSE()
        target_link_libraries(${testName} ${Boost_LIBRARIES} ${LIBS})
    ENDIF()

    #Finally add it to test execution -
    #Notice the WORKING_DIRECTORY and COMMAND
    add_test(NAME ${testName} COMMAND ${testName} )
ELSE()
    MESSAGE ("Did not find Boost unit_test_framework. Will compile unit tests.")
ENDIF()
