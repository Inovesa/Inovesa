/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2013-2016: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2016: Karlsruhe Institute of Technology                 *
 *                                                                            *
 * This file is part of Inovesa.                                              *
 * Inovesa is free software: you can redistribute it and/or modify            *
 * it under the terms of the GNU General Public License as published by       *
 * the Free Software Foundation, either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * Inovesa is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU General Public License for more details.                               *
 *                                                                            *
 * You should have received a copy of the GNU General Public License          *
 * along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           *
 ******************************************************************************/

#ifndef PHASESPACE_HPP
#define PHASESPACE_HPP

#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <GL/glew.h>
#ifndef __APPLE__
#include <GL/gl.h>
#else // non-Apple
#include <OpenGL/gl.h>
#endif // non-Apple
#include <list>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace vfps {
        class PhaseSpace; // forward declaration
}

#include "CL/OpenCLHandler.hpp"
#include "defines.hpp"
#include "Ruler.hpp"

namespace vfps
{

class PhaseSpace
{
public:
    struct Position {
        meshaxis_t x;
        meshaxis_t y;
    };

public:
    enum class IntegralType : uint_fast8_t {
        sum,simpson
    };

public:
    PhaseSpace() = delete;

    PhaseSpace(std::array<Ruler<meshaxis_t>,2> axis,
               const double bunch_charge, const double bunch_current,
               const double zoom=1, meshdata_t *data = nullptr);

    PhaseSpace(Ruler<meshaxis_t> axis1, Ruler<meshaxis_t> axis2,
               const double bunch_charge, const double bunch_current,
               const double zoom=1, meshdata_t *data = nullptr);

    PhaseSpace(meshindex_t ps_size,
               meshaxis_t xmin, meshaxis_t xmax,
               meshaxis_t ymin, meshaxis_t ymax,
               const double bunch_charge, const double bunch_current,
               double xscale=0, double yscale=0,
               const double zoom=1, meshdata_t *data = nullptr);

    PhaseSpace(const PhaseSpace& other);

    ~PhaseSpace() noexcept;

     /**
      * @brief getData gives direct access to held data
      *
      * @return pointer to array holding size<0>()*size<1>() data points
      */
    inline meshdata_t* getData() const
    { return _data1D; }

    inline meshaxis_t getDelta(const uint_fast8_t x) const
    { return _axis[x].delta(); }

    inline meshaxis_t getMax(const uint_fast8_t x) const
    { return _axis[x].max(); }

    inline meshaxis_t getMin(const uint_fast8_t x) const
    { return _axis[x].min(); }

    inline double getScale(const uint_fast8_t x) const
    { return _axis[x].scale(); }

    /**
     * @brief getAxis
     * @param x which axis? (0 -> x or 1 -> y)
     * @return reference to the Axis describing mesh in x direction
     */
    inline const Ruler<meshaxis_t>* getAxis(const uint_fast8_t x) const
    { return &(_axis[x]); }

    /**
     * @brief average
     * @param axis which axis? (0 -> x or 1 -> y)
     * @return projection to axis
     *
     * relies on an up-t date _projection[axis]
     */
    meshdata_t average(const uint_fast8_t axis);

    /**
     * @brief integral
     * @return integrated phase space volume
     *
     * relies on an up-to-date _projection[0]
     */
    integral_t integral();

    const integral_t& getIntegral() const
    { return _integral; }

    /**
     * @brief variance
     * @param axis which axis? (0 -> x or 1 -> y)
     * @return bunch length / energy spread
     *
     * relies on an up-to-date _projection[axis]
     */
    meshdata_t variance(const uint_fast8_t axis);

    inline meshdata_t getBunchLength() const
        { return std::sqrt(getMoment(0,1)); }

    inline meshdata_t getEnergySpread() const
        { return std::sqrt(getMoment(1,1)); }

    inline meshdata_t getMoment(const uint_fast8_t x,const uint_fast16_t m) const
        { return _moment[x][m]; }

    inline projection_t* getProjection(const uint_fast8_t x) const
        { return _projection[x]; }

    /**
     * @brief updateXProjection
     *
     * @todo make ready for arbitrary meshdata_t (currently hardcoded to float)
     */
    void updateXProjection();

    void updateYProjection();

    /**
     * @brief normalize
     * @return integral before normalization
     *
     * @todo: Use OpenCL
     *
     * normalize() does neither recompute the integral nor sets it to 1
     */
    integral_t normalize();

    inline meshdata_t* operator[](const meshindex_t i) const
    { return _data[i]; }

    PhaseSpace& operator=(PhaseSpace other);

    inline size_t nMeshCells() const
    { return _axis[0].steps()*_axis[1].steps(); }

    inline size_t nMeshCells(const uint_fast8_t x) const
    { return _axis[x].steps(); }

    inline meshaxis_t size(const uint_fast8_t x) const
    { return _axis[x].size(); }

    inline meshaxis_t x(const uint_fast8_t axis, const size_t n) const
        { return _axis[axis][n]; }

    /**
     * @brief swap
     * @param other
     */
    friend void swap(PhaseSpace& first, PhaseSpace& second) noexcept;

    #ifdef INOVESA_USE_CL
    void syncCLMem(clCopyDirection dir);
    #endif

protected:
    const std::array<Ruler<meshaxis_t>,2> _axis;

public:
    /**
     * @brief _charge conversion factor _integral -> bunch charge in C
     */
    const double charge;

    /**
     * @brief _charge conversion factor _integral -> bunch current in A
     */
    const double current;

protected:
    /**
     * @brief _integral as we work in normalitzed units, this should be 1
     */
    integral_t _integral;

    const std::array<projection_t*,2> _projection;

    const uint32_t _nmeshcellsX;

    const uint32_t _nmeshcellsY;

    const size_t _nmeshcells;

    const IntegralType _integraltype;

    meshdata_t** _data;

    meshdata_t* _data1D;

    /**
     * @brief _moment: holds the moments for distributions
     *            in both axis in mesh coordinates
     *
     * 0: average
     * 1: variance
     * 2: skewness
     * 3: kurtosis
     */
    std::array<std::array<meshdata_t,4>,2> _moment;

    meshdata_t* _ws;

#ifdef INOVESA_USE_CL
public:
    cl::Buffer data_buf;

    cl::Buffer projectionX_buf;

    cl::Buffer integral_buf;

private:
    cl::Program _clProgProjX;

    cl::Kernel  _clKernProjX;

    cl::Program _clProgIntegral;

    cl::Kernel  _clKernIntegral;

    cl::Buffer  ws_buf;

    static std::string cl_code_integral;

    static std::string cl_code_projection_x;
#endif // INOVESA_USE_CL

public:
    void createFromProjections();

    /**
     * @brief gaus calculates gaussian distribution
     * @param axis
     * @param zoom
     */
    void gaus(const uint_fast8_t axis, const double zoom);
};

/**
 * @brief swap
 * @param first
 * @param second
 *
 * @todo adjust to also swap cl::Buffer
 */
void swap(PhaseSpace& first, PhaseSpace& second) noexcept;

}

#endif // PHASESPACE_HPP
