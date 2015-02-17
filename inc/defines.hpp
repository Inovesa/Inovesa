#ifndef DEFINES_HPP
#define DEFINES_HPP

#include "fixed_point.h"

#define INOVESA_USE_GUI
#define INOVESA_USE_CL

namespace vfps
{

constexpr double rotations = 1;
constexpr unsigned int patterndim_x = 512;
constexpr unsigned int patterndim_y = 4048;
constexpr unsigned int pattern_margin = 128;

enum class pattern {
	square, gaus, half, quarters
};
constexpr pattern ptrntype = pattern::quarters;

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

typedef float meshaxis_t;
typedef float meshdata_t;
typedef float interpol_t;
typedef float integral_t;
}

#endif // DEFINES_HPP
