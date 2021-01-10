// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include <memory>
#include <sstream>

#include "defines.hpp"
#include "IO/Display.hpp"
#include "PS/PhaseSpace.hpp"

namespace vfps
{

class SourceMap
{
public:
    enum InterpolationType : uint_fast8_t {
        none = 1,
        linear = 2,
        quadratic = 3,
        cubic = 4
    };

    enum class RotationCoordinates : uint_fast8_t {
        phys_pq = 0,
        mesh = 1, // rotate on mesh
        norm_0_1 = 2, // normalized space between 0 and 1
        norm_pm1 = 3 // normalized space between -1 and +1
    };

protected:
    typedef struct {
        meshindex_t index;
        interpol_t weight;
    } hi;

public:
    /**
     * @brief SourceMap
     * @param in
     * @param out
     * @param xsize
     * @param ysize
     * @param memsize number of hi (needed for memory optimization)
     * @param interpoints
     * @param intertype number of points used for interpolation
     */
    SourceMap( std::shared_ptr<PhaseSpace> in, std::shared_ptr<PhaseSpace> out
             , meshindex_t xsize, meshindex_t ysize, size_t memsize
             , uint_fast8_t interpoints, uint_fast8_t intertype
             , oclhptr_t oclh
             );

    /**
     * @brief SourceMap
     * @param in
     * @param out
     * @param xsize
     * @param ysize
     * @param interpoints number of points used for interpolation
     */
    SourceMap( std::shared_ptr<PhaseSpace> in
             , std::shared_ptr<PhaseSpace> out
             , size_t xsize, size_t ysize
             , uint_fast8_t interpoints, uint_fast8_t intertype
             , oclhptr_t oclh
             );

    /**
     * @brief ~SourceMap deletes OpenCL events (if INOVESA_ENABLE_CLPROFILING == 1)
     */
    virtual ~SourceMap() noexcept;

    /**
     * @brief apply the SM
     */
    virtual void apply();

    /**
     * @brief apply
     * @param pos
     */
    virtual void applyTo(PhaseSpace::Position& pos) const=0;

    void applyToAll(std::vector<PhaseSpace::Position>& particles);

protected:
    /**
     * @brief _ip holds the total number of points used for interpolation
     */
    #if INOVESA_USE_OPENCL == 1
    const cl_uint _ip;
    #else
    const uint_fast8_t _ip;
    #endif // INOVESA_USE_OPENCL

    /**
     * @brief _ip holds the per dimension number
     *            of points used for interpolation
     */
    const uint_fast8_t _it;

    /**
     * @brief _hinfo
     */
    std::vector<hi> _hinfo;

    /**
     * @brief _xsize horizontal size of the SourceMap (in grid points)
     */
    const meshindex_t _xsize;

    /**
     * @brief _ysize vertical size of the SourceMap (in grid points)
     */
    const meshindex_t _ysize;

    #if INOVESA_USE_OPENCL == 1
    /**
     * @brief _hi_buf buffer for source information
     */
    cl::Buffer _sm_buf;

    /**
     * @brief applySM
     */
    cl::Kernel applySM;

    #if INOVESA_ENABLE_CLPROFILING == 1
    std::unique_ptr<cl::vector<cl::Event*>> applySMEvents;

    std::unique_ptr<cl::vector<cl::Event*>> syncSMEvents;
    #endif // INOVESA_ENABLE_CLPROFILING

    std::string _cl_code;

    cl::Program _cl_prog;

    #endif // INOVESA_USE_OPENCL

    std::array<meshRuler_ptr,2> _axis;

    std::shared_ptr<PhaseSpace> _in;
    std::shared_ptr<PhaseSpace> _out;

    oclhptr_t _oclh;

    #if INOVESA_USE_OPENCL == 1
    /**
     * @brief genCode4SM1D generates OpenCL code for a generic source map
     */
    void genCode4SM1D();
    #endif // INOVESA_USE_OPENCL

    /**
     * @brief calcCoefficiants
     * @param ic array to store interpolation coefficiants
     * @param f distance from lower mesh point
     * @param it number of interpolation coefficiants (size of ic)
     */
    static void calcCoefficiants(interpol_t* ic, const interpol_t f,
                          const uint_fast8_t it);

    static void notClampedMessage();

    #if INOVESA_ENABLE_CLPROFILING == 1
    void saveTimings(std::string mapname);
    #endif // INOVESA_ENABLE_CLPROFILING
};

}
