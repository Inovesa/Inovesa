// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

#include "SM/WakeKickMap.hpp"

namespace vfps
{

/**
 * @brief The WakeFunctionMap class
 */
class WakeFunctionMap : public WakeKickMap
{
public:
    WakeFunctionMap() = delete;

    WakeFunctionMap( std::shared_ptr<PhaseSpace> in
                   , std::shared_ptr<PhaseSpace> out
                   , const meshindex_t xsize, const meshindex_t ysize
                   , const std::string& fname
                   , const double sigmaE, const double E0
                   , const double Ib, const double dt
                   , const InterpolationType it
                   , const bool interpol_clamp
                   , oclhptr_t oclh
                   );

    WakeFunctionMap( std::shared_ptr<PhaseSpace> in
                   , std::shared_ptr<PhaseSpace> out
                   , const meshindex_t xsize, const meshindex_t ysize
                   , const ElectricField* csr
                   , const InterpolationType it
                   , const bool interpol_clamp
                   , oclhptr_t oclh
                   );

    ~WakeFunctionMap() noexcept;

public:
    inline const meshaxis_t* getWakeFunction() const
        { return _wakefunction; }

    /**
     * @brief update overrides WakeKickMap
     *
     * @todo currently uses phasespace in CPU Ram
     */
    void update() override;

private:
    /**
     * @brief _wakeFromFile reads in a file and scales wake to internal units
     * @param fname file name to read wake from
     */
    void _wakeFromFile(const std::string fname, const double scaling);


private:
    /**
     * @brief WakeFunctionMap
     * @param in
     * @param out
     * @param xsize
     * @param ysize
     * @param it
     * @param interpol_clamp
     * @param oclh
     *
     * @todo currently broken when used with CL/GL sharing
     */
    WakeFunctionMap( std::shared_ptr<PhaseSpace> in
                   , std::shared_ptr<PhaseSpace> out
                   , const meshindex_t xsize, const meshindex_t ysize
                   , const InterpolationType it, bool interpol_clamp
                   , oclhptr_t oclh
                   );


    const Ruler<meshaxis_t> _xaxis;

    /**
     * @brief _wakefunktion (normalized single particle) wake,
     *          sampled at 2*xsize positions [-xsize:+xsize]
     */
    meshaxis_t* const _wakefunction;

    const size_t _wakesize;
};

} // namespace vfps

