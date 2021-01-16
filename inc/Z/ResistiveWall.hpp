// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include "Z/Impedance.hpp"

namespace vfps
{

/**
 * @brief The ResistiveWall class models the resistive wall impedance
 */
class ResistiveWall : public Impedance
{
public:
    /**
     * @brief ResistiveWall
     * @param n
     * @param f0
     * @param f_max
     * @param L: Length of the beam pipe
     * @param s: conductivity [S/m]
     * @param xi
     * @param b
     */
    ResistiveWall( const size_t n
                 , const frequency_t f0
                 , const frequency_t f_max
                 , const double L
                 , const double s
                 , const double xi
                 , const double b
                 , oclhptr_t oclh = nullptr
                 );

private:
    static std::vector<vfps::impedance_t>
    __calcImpedance(const size_t n,
                    const frequency_t f0,
                    const frequency_t f_max,
                    const double L,
                    const double s,
                    const double xi,
                    const double b);
};

} // namespace vfps

