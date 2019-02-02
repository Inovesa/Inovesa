// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "CL/OpenCLHandler.hpp"
#include "defines.hpp"
#include "PS/Ruler.hpp"

namespace vfps
{

class Impedance
{
public:
    Impedance() = delete;

    /**
     * @brief Impedance copy constructor
     * @param other
     */
    Impedance(const Impedance &other);

    /**
     * @brief Impedance basic constructor that initializes everything
     * @param axis
     * @param z
     */
    Impedance(Ruler<frequency_t> &&axis
             , const std::vector<impedance_t> &z
             , oclhptr_t oclh = nullptr
             );


    /**
     * @brief Impedance
     * @param z
     * @param f_max
     *
     * Note that we will use this for DFT,
     * so n>z.size()/2 is defined to be equivalent to n<0.
     */
    Impedance( const std::vector<impedance_t>& z
             , const frequency_t f_max
             , oclhptr_t oclh = nullptr
             );


    Impedance( const size_t nfreqs
             , const frequency_t f_max
             , oclhptr_t oclh = nullptr
             );


    /**
     * @brief Impedance
     * @param name of datafile in the format "n Re(Z) Im(Z)",
     *        where n=f/f_rev is the revolution harmonic
     */
    Impedance( std::string datafile, double f_max, oclhptr_t oclh = nullptr);

    inline const impedance_t* data() const
        { return _data.data(); }

    inline const std::vector<impedance_t>& impedance() const
        { return _data; }

    inline const impedance_t operator[](size_t n) const
        { return _data[n]; }

    inline size_t nFreqs() const
        { return _nfreqs; }

    inline size_t size() const
        { return _data.size(); }

    inline const Ruler<frequency_t>* getRuler() const
        { return &_axis; }

    #if INOVESA_USE_OPENCL == 1
    cl::Buffer data_buf;
    #endif // INOVESA_USE_OPENCL

    static constexpr double factor4Ohms = 1;

    /// vacuum impedance (in Ohms)
    static constexpr double Z0 = 376.730313461;

public:
    /**
     * @brief operator +=
     * @param rhs
     * @return
     *
     * assumes size() equals rhs.size()
     */
    Impedance& operator+=(const Impedance& rhs);

private:
    size_t _nfreqs;

    const Ruler<frequency_t> _axis;

protected:
    std::vector<impedance_t> _data;

    void syncCLMem();

    oclhptr_t _oclh;

private:
    static std::vector<impedance_t> readData(std::string fname);

public:
    /**
     * @brief upper_power_of_two
     * @param v
     * @return
     *
     * @todo In Impedance for historical reasons.
     * It might be useful to move it somewhere else,
     * as it does now relate to more than just the impedance.
     *
     * see http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
     */
    static uint64_t upper_power_of_two(uint64_t v);
};

inline Impedance operator+(const Impedance& lhs, const Impedance& rhs)
{
    return Impedance(lhs) += rhs;
}

} // namespace vfps

