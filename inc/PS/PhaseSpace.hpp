// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include <algorithm>
#include <array>
#include <boost/multi_array.hpp>
#include <cfloat>
#include <cmath>
#include <fstream>
#if INOVESA_USE_OPENGL == 1
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
    PhaseSpace() = delete;

    /**
     * PhaseSpace initalizing constructor
     */
    PhaseSpace(std::array<meshRuler_ptr,2> axis
              , oclhptr_t oclh
              , const double beam_charge
              , const double beam_current
              , const std::vector<integral_t>& filling={1}
              , const double zoom=1
              , const meshdata_t *data = nullptr
              );

    /**
     * PhaseSpace initalizing constructor
     */
    PhaseSpace(meshRuler_ptr axis0
              , meshRuler_ptr axis1
              , oclhptr_t oclh
              , const double beam_charge
              , const double beam_current
              , const std::vector<integral_t>& filling={1}
              , const double zoom=1
              , const meshdata_t* data = nullptr
              );

    /**
     * @brief PhaseSpace initalizing constructor
     */
    PhaseSpace( meshaxis_t qmin
              , meshaxis_t qmax
              , double qscale
              , meshaxis_t pmin
              , meshaxis_t pmax
              , double pscale
              , oclhptr_t oclh
              , const double beam_charge
              , const double beam_current
              , const std::vector<integral_t>& filling={1}
              , const double zoom=1
              , const meshdata_t* data = nullptr
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
    inline const meshdata_t* getData() const
    { return _data.data(); }

    inline meshdata_t* getData()
    { return _data.data(); }

    inline auto operator [] (const unsigned int i)
    { return _data[i]; }

    inline meshaxis_t getDelta(const uint_fast8_t x) const
    { return _axis[x]->delta(); }

    inline meshaxis_t getMax(const uint_fast8_t x) const
    { return _axis[x]->max(); }

    inline meshaxis_t getMin(const uint_fast8_t x) const
    { return _axis[x]->min(); }

    inline const auto& getScale(const uint_fast8_t x) const
    { return _axis[x]->scale(); }

    inline auto getScale(const uint_fast8_t x, std::string unit) const
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

    inline const std::vector<integral_t> getBunchPopulation() const
        { return _filling; }

    inline auto getSetBunchPopulation() const
        { return _filling_set; }

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

    /**
     * @brief getMoment
     * @param x axis
     * @param m m-th moment
     */
    inline auto getMoment( const uint_fast8_t x
                         , const uint_fast8_t m) const
        { return _moment[x][m]; }

    inline auto getBunchLength() const
        { return _rms[0]; }

    inline auto getEnergySpread() const
        { return _rms[1]; }


    inline auto getProjection(const uint_fast8_t x) const
        { return _projection[x]; }

    /**
     * @brief updateXProjection updates longitudinal bunch profiles
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
    inline const std::vector<integral_t>& integrateAndNormalize() {
        integrate();
        normalize();
        return _filling;
    }

    const std::vector<integral_t>& normalize();

    /**
     * @brief operator = unifying assignment operator
     * @param other
     * @return
     */
    PhaseSpace& operator =(PhaseSpace other);

    /**
     * @brief nBunches number of RF buckets in simulation
     * @return
     */
    inline meshindex_t nBunches() const
    { return _nbunches; }

    /**
     * @brief nMeshCells total number of mesh cells
     * @return
     */
    inline meshindex_t nMeshCells() const
    { return _axis[0]->steps()*_axis[1]->steps(); }

    /**
     * @brief nMeshCells number of mes cells in direction
     * @param x direction (0: x, 1: y)
     * @return
     */
    inline meshindex_t nMeshCells(const uint_fast8_t x) const
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
     * @param other
     *
     * @todo adjust to also swap cl::Buffer and other elements
     */
    void swap(PhaseSpace& other) noexcept;

    #if INOVESA_USE_OPENCL == 1
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

public:
    /**
     * @brief nx reference to _nmeshcellsX
     */
    static const meshindex_t& nx;

    /**
     * @brief ny reference to _nmeshcellsY
     */
    static const meshindex_t& ny;

    /**
     * @brief nb reference to _nbunches
     */
    static const meshindex_t& nb;

    /**
     * @brief nxy reference to _nmeshcells
     */
    static const meshindex_t& nxy;

    /**
     * @brief nxy reference to _totalmeshcells
     */
    static const meshindex_t& nxyb;


    #if INOVESA_ALLOW_PS_RESET == 1
    /**
     * @brief resetSize relevant for unit tests
     */
    inline static void resetSize()
        { _firstinit = true; }


    static void resetSize( const meshindex_t x,
                           const meshindex_t b)
    {
        resetSize();
        setSize(x,b);
    }

    #endif

    /**
     * @brief setSize one-time setter for sizes
     * @param x number of grid cells in x (and y) directions
     * @param b number of bunches
     *
     * As all grids have to have the same size, it is set globally.
     *
     * Support for non-quadratic grids is considered for a future release.
     */
    static void setSize(const meshindex_t x,
                        const meshindex_t b);

protected:
    static bool _firstinit;

    /**
     * @brief _nmeshcellsX number of cells for position axis
     */
    static meshindex_t _nmeshcellsX;

    /**
     * @brief _nmeshcellsY number of cells for energy axis
     */
    static meshindex_t _nmeshcellsY;

    /**
     * @brief _nbunches number of bunches (phase spaces)
     */
    static meshindex_t _nbunches;

    /**
     * @brief _nmeshcells number of cells for one phase space
     */
    static meshindex_t _nmeshcells;

    /**
     * @brief _totalmeshcells accumuated number of cells for all phase spaces
     */
    static meshindex_t _totalmeshcells;

protected:
    /**
     * @brief _fillingpattern: normalized bunch charges as they should be
     *
     * as we work in normalitzed units, the sum should be 1, e.g. {0.25,0.75},
     * empty buckets are omited
     */
    const std::vector<integral_t> _filling_set;

    /**
     * @brief _filling: normalized bunch charges as they are
     *
     * ideally, this is the same as _filling_set
     */
    std::vector<integral_t> _filling;

   /**
    * @brief _integral
    *
    * as we work in normalitzed units, this should be 1
    */
   integral_t _integral;

    /**
     * @brief _projection dimensions are orientation, bunch, x/y grid cell
     */
    boost::multi_array<projection_t,3> _projection;

    /**
     * @brief _data dimensions are: bunch, x coordinate, y coordinate
     */
    boost::multi_array<meshdata_t,3> _data;

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
    boost::multi_array<meshaxis_t,3> _moment;


    /**
     * @brief _RMS
     *
     * 0: x-Axis
     * 1: y-axis
     *
     * n: bunch number
     */
    boost::multi_array<meshaxis_t,2> _rms;

    /**
     * @brief _ws weights for Simpson integration
     *
     * assumes nx == ny
     */
    const std::vector<meshdata_t> _ws;

private:
    oclhptr_t _oclh;

public:
    #if INOVESA_USE_OPENGL == 1
    /**
     * @brief buffer to show projection in OpenGL view
     *
     * projectionX_glbuf will be allocated by the constructor
     * if CL/GL sharing is enabled. Else, it will just hold
     * 0 and allocation is left to the Inovesa displaying system.
     */
    vfps::clgluint projectionX_glbuf;
    #endif // INOVESA_USE_OPENGL

    #if INOVESA_USE_OPENCL == 1
    cl::Buffer data_buf;

    cl::Buffer projectionX_clbuf;

    cl::Buffer filling_buf;

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
    #endif // INOVESA_USE_OPENCL == 1

private:
    void createFromProjections();

    /**
     * @brief gaus calculates gaussian distribution
     * @param axis
     * @param zoom
     */
    void gaus(const uint_fast8_t axis,
              const meshindex_t bunch,
              const double zoom);

    /**
     * @brief simpsonWeights helper function to allow for const _ws
     * @return
     */
    const std::vector<meshdata_t> simpsonWeights();
};

void swap(PhaseSpace& first, PhaseSpace& second) noexcept;

} // namespace vfps

namespace std {
/**
 * @brief specialization of std::swap to use the above vfps::swap
 */
template<>
inline void swap(vfps::PhaseSpace& first, vfps::PhaseSpace& second)
    {vfps::swap(first, second); }

} // namespace std
