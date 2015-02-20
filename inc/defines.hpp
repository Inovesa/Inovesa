#ifndef DEFINES_HPP
#define DEFINES_HPP

#include "fixed_point.h"

#define INOVESA_USE_GUI
#define INOVESA_USE_CL

/**
 * possible choices are:
 * 1: sum
 * 2: simpson
 */
#define INTEGRAL_TYPE 2

/**
  * possible choices are:
  * 1: no interpolation
  * 2: linear interpolation
  * 3: quadratic interpolation
  * 4: cubic interpolation
  */
#define INTERPOL_TYPE 4

/**
 * possible choices are:
 * 1: rotate on mesh
 * 2: normalized space between 0 and 1
 * 3: normalized space between -1 and +1
 */
#define ROTATION_TYPE 3

/**
 * possible choices are:
 * 0: no test pattern
 * 1: a square
 * 2: a 2D gaussian
 * 3: a rectengle (half)
 * 4: quarters with different patterns
 */
#define TEST_PATTERN 4

namespace vfps
{

constexpr double rotations = 1;
constexpr unsigned int patterndim_x = 512;
constexpr unsigned int patterndim_y = 4048;
constexpr unsigned int pattern_margin = 128;

constexpr unsigned int steps = 4000;

/**
 * @brief ps_xsize horizontal size of the phase space (in mesh points)
 */
constexpr unsigned int ps_xsize = 512;

/**
 * @brief ps_ysize vertical size of the phase space (in mesh points)
 */
constexpr unsigned int ps_ysize = 512;

typedef fpml::fixed_point<int32_t,2,29> fixp32;
typedef fpml::fixed_point<int64_t,34,29> fixp64;

typedef float meshaxis_t;
typedef fixp32 meshdata_t;
typedef fixp32 interpol_t;
typedef fixp64 integral_t;
}

#endif // DEFINES_HPP
