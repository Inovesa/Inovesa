/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   *
 * Copyright (c) 2014-2015: Patrik Schönfeldt                                 *
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

#ifndef IMPEDANCE_HPP
#define IMPEDANCE_HPP

#include <fstream>
#include <string>
#include <vector>

#include "defines.hpp"
#include "Ruler.hpp"

namespace vfps
{

class Impedance
{
public:
	enum class ImpedanceModel : uint_fast8_t {
		FreeSpace = 0
	};

public:
	/**
	 * @brief Impedance
	 * @param model
	 * @param n compute (at least) to wavenumber n=f/f_rev,
	 *        where f_rev is the revolution frequency
	 * @param roundup round to next 2^N
	 */
	Impedance(ImpedanceModel model, size_t n, bool roundup=true);

	/**
	 * @brief Impedance
	 * @param name of datafile in the format "n Re(Z) Im(Z)",
	 *        where n=f/f_rev is the wavenumber
	 */
	Impedance(std::string datafile);

    inline const impedance_t* data() const
		{ return _data.data(); }

    inline const std::vector<impedance_t>& impedance() const
        { return _data; }

	inline const impedance_t operator[](size_t n) const
		{ return _data[n]; }

	inline size_t maxN() const
		{ return _data.size(); }

private:
	std::vector<impedance_t> _data;

	size_t _nmax;

	/**
	 * @brief upper_power_of_two
	 * @param v
	 * @return
	 *
	 * see http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
	 */
	uint64_t upper_power_of_two(uint64_t v);


};

} // namespace vfps

#endif // IMPEDANCE_HPP
