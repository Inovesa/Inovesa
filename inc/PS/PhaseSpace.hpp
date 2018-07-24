/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2013-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
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
#ifdef INOVESA_USE_OPENGL
#include <GL/glew.h>
#ifndef __APPLE__
#include <GL/gl.h>
#else // non-Apple
#include <OpenGL/gl.h>
#endif // non-Apple
#endif // INOVESA_USE_OPENGL
#include <list>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <vector>

namespace vfps {
        class PhaseSpace; // forward declaration
}

#include "Array.h"

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
    enum class IntegralMethod : uint_fast8_t {
        sum,simpson
    };

public:
    PhaseSpace() = delete;

    PhaseSpace( std::array<meshRuler_ptr,2> axis
              , oclhptr_t oclh
              , const double bunch_charge
              , const double bunch_current
              , const uint32_t nbunches=1
              , const double zoom=1
              , meshdata_t* data = nullptr
              );

    PhaseSpace( meshRuler_ptr axis0
              , meshRuler_ptr axis1
              , oclhptr_t oclh
              , const double bunch_charge
              , const double bunch_current
              , const uint32_t nbunches=1
              , const double zoom=1
              , meshdata_t* data = nullptr
              );

    PhaseSpace( meshindex_t ps_size
              , meshaxis_t qmin
              , meshaxis_t qmax
              , double qscale
              , meshaxis_t pmin
              , meshaxis_t pmax
              , double pscale
              , oclhptr_t oclh
              , const double bunch_charge
              , const double bunch_current
              , const uint32_t nbunches=1
              , const double zoom=1
              , meshdata_t *data = nullptr
              );

    /**
     * @brief PhaseSpace copy constructor
     */
    PhaseSpace(const PhaseSpace& other);

    ~PhaseSpace() noexcept;

     /**
      * @brief getData gives direct access to held data
      *
      * @return pointer to array holding size<0>()*size<1>() data points
      */
    inline meshdata_t* getData() const
    { return _data(); }

    inline auto operator [] (const unsigned int i)
    { return _data[i]; }

    inline meshaxis_t getDelta(const uint_fast8_t x) const
    { return _axis[x]->delta(); }

    inline meshaxis_t getMax(const uint_fast8_t x) const
    { return _axis[x]->max(); }

    inline meshaxis_t getMin(const uint_fast8_t x) const
    { return _axis[x]->min(); }

    inline std::map<std::string,double> getScale(const uint_fast8_t x) const
    { return _axis[x]->scale(); }

    inline double getScale(const uint_fast8_t x, std::string unit) const
    { return _axis[x]->scale(unit); }

    /**
     * @brief getAxis
     * @param x which axis? (0 -> x or 1 -> y)
     * @return reference to the Axis describing mesh in x direction
     */
    inline const meshRuler_ptr getAxis(const uint_fast8_t x) const
    { return _axis[x]; }

    /**
     * @brief average
     * @param axis which axis? (0 -> x or 1 -> y)
     *
     * relies on an up-t date _projection[axis]
     */
    void average(const uint_fast8_t axis);

    /**
     * @brief integral
     * @return integrated phase space volume
     *
     * relies on an up-to-date _projection[0]
     */
    void integrate();

    inline const Array::array1<integral_t> getBunchPopulation() const
        { return _bunchpopulation; }

    inline integral_t getIntegral() const
        { return _integral; }

    /**
     * @brief variance
     * @param axis which axis? (0 -> x or 1 -> y)
     * @return bunch length / energy spread
     *
     * relies on an up-to-date _projection[axis]
     */
    void variance(const uint_fast8_t axis);

    inline auto getMoment( const uint_fast8_t x
                         , const uint_fast8_t m) const
        { return _moment[x][m]; }

    inline const meshaxis_t* getBunchLength() const
        { return _rms[0]; }

    inline const meshaxis_t* getEnergySpread() const
        { return _rms[1]; }


    inline const projection_t* getProjection(const uint_fast8_t x) const
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
    Array::array1<integral_t> normalize();

    PhaseSpace& operator=(PhaseSpace other);

    /**
     * @brief nBunches number of RF buckets in simulation
     * @return
     */
    inline uint32_t nBunches() const
    { return _nbunches; }

    /**
     * @brief nMeshCells total number of mesh cells
     * @return
     */
    inline size_t nMeshCells() const
    { return _axis[0]->steps()*_axis[1]->steps(); }

    /**
     * @brief nMeshCells number of mes cells in direction
     * @param x direction (0: x, 1: y)
     * @return
     */
    inline uint32_t nMeshCells(const uint_fast8_t x) const
    { return _axis[x]->steps(); }

    /**
     * @brief length in reduced coordinates
     * @param x
     * @return
     */
    inline meshaxis_t length(const uint_fast8_t x) const
    { return _axis[x]->length(); }

    /**
     * @brief x grid coordinate of normalized position
     * @param q
     * @return
     */
    inline meshaxis_t x(const meshaxis_t q) const
        { return std::min( std::max( 0.0f,(q-_axis[0]->min())/_axis[0]->delta() )
                         , _nmeshcellsX-1.0f) ; }

    /**
     * @brief y grid coordinate of normalized energy
     * @param p
     * @return
     */
    inline meshaxis_t y(const meshaxis_t p) const
        { return std::min( std::max( 0.0f,(p-_axis[1]->min())/_axis[1]->delta() )
                         , _nmeshcellsY-1.0f) ; }

    /**
     * @brief q normalized position coordinate of grid point x
     * @param x
     * @return
     */
    inline meshaxis_t q(const meshindex_t x) const
        { return _qp(0,x); }

    /**
     * @brief p normalized energy coordinate of grid point y
     * @param y
     * @return
     */
    inline meshaxis_t p(const meshindex_t y) const
        { return _qp(1,y); }

private:
    inline meshaxis_t _qp(const uint_fast8_t axis, const meshindex_t n) const
        { return _axis[axis]->at(n); }

public:
    /**
     * @brief swap
     * @param first
     * @param second
     */
    friend void swap(PhaseSpace& first, PhaseSpace& second) noexcept;

    #ifdef INOVESA_USE_OPENCL
    void syncCLMem(OCLH::clCopyDirection dir, cl::Event* evt = nullptr);
    #endif // INOVESA_USE_OPENCL

protected:
    const std::array<meshRuler_ptr,2> _axis;

public:
    /**
     * @brief charge conversion factor _integral -> bunch charge in C
     */
    const double charge;

    /**
     * @brief current conversion factor _integral -> bunch current in A
     */
    const double current;

protected:
    const uint32_t _nmeshcellsX;

    const uint32_t _nmeshcellsY;

    const uint32_t _nbunches;

    const size_t _nmeshcells;

    const IntegralMethod _integralmethod;

    /**
     * @brief _bunchpopulation as we work in normalitzed units,
     * the sum should be 1
     */
    Array::array1<integral_t> _bunchpopulation;

    /**
     * @brief _integral as we work in normalitzed units, this sum should be 1
     */
    integral_t _integral;

    /**
     * @brief _projection dimensions are orientation, bunch, x/y grid cell
     */
    Array::array3<projection_t> _projection;

    /**
     * @brief _data dimensions are: bunch, x coordinate, y coordinate
     */
    Array::array3<meshdata_t> _data;

    /**
     * @brief _moment: holds the moments for distributions
     *            in both axis in mesh coordinates
     *
     * 0: x-Axis
     * 1: y-axis
     *
     * 0: average
     * 1: variance
     * 2: skewness
     * 3: kurtosis
     *
     * n: bunch number
     */
    Array::array3<meshaxis_t> _moment;


    /**
     * @brief _RMS
     *
     * 0: x-Axis
     * 1: y-axis
     *
     * n: bunch number
     */
    Array::array2<meshaxis_t> _rms;

    /**
     * @brief _ws weights for Simpson integration
     */
    const Array::array1<meshdata_t> _ws;

private:
    oclhptr_t _oclh;

public:
    #ifdef INOVESA_USE_OPENGL
    /**
     * @brief projectionX_glbuf will only be allocated during CL/GL sharing
     */
    cl_GLuint projectionX_glbuf;
    #endif // INOVESA_USE_OPENGL
    #ifdef INOVESA_USE_OPENCL

    cl::Buffer data_buf;

    cl::Buffer projectionX_clbuf;

    cl::Buffer bunchpop_buf;

private:
    cl::Program _clProgProjX;

    cl::Kernel  _clKernProjX;

    cl::Program _clProgIntegral;

    cl::Kernel  _clKernIntegral;

    std::unique_ptr<cl::vector<cl::Event*>> xProjEvents;

    std::unique_ptr<cl::vector<cl::Event*>> integEvents;

    std::unique_ptr<cl::vector<cl::Event*>> syncPSEvents;

    cl::Buffer  ws_buf;

    static std::string cl_code_integral;

    static std::string cl_code_projection_x;
#endif // INOVESA_USE_OPENCL

private:
    void createFromProjections();

    /**
     * @brief gaus calculates gaussian distribution
     * @param axis
     * @param zoom
     */
    void gaus(const uint_fast8_t axis, const uint32_t bunch, const double zoom);

    /**
     * @brief simpsonWeights helper function to allow for const _ws
     * @return
     */
    const Array::array1<meshdata_t> simpsonWeights();
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
