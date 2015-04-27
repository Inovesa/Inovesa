/******************************************************************************/
/* Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   */
/* Copyright (c) 2014-2015: Patrik Schönfeldt                                 */
/*                                                                            */
/* This file is part of Inovesa.                                              */
/* Inovesa is free software: you can redistribute it and/or modify            */
/* it under the terms of the GNU General Public License as published by       */
/* the Free Software Foundation, either version 3 of the License, or          */
/* (at your option) any later version.                                        */
/*                                                                            */
/* Inovesa is distributed in the hope that it will be useful,                 */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           */
/******************************************************************************/

#ifndef KICKMAP_HPP
#define KICKMAP_HPP

#include <array>

#include "defines.hpp"
#include "HeritageMap.hpp"

namespace vfps
{

/**
 * @brief The KickMap class offers an option for one-dimensional kicks
 *
 * @todo implement option to have 1D interpolation with HeritageMap
 */
class WakeKickMap : public HeritageMap
{
public:
	/**
	 * @brief WakeKickMap
	 * @param in
	 * @param out
	 * @param xsize
	 * @param ysize
	 * @param wake has to be normalized in a way,
	 *	      that plain multiplication with density gives the force
	 */
	WakeKickMap(PhaseSpace* in, PhaseSpace* out,
			const unsigned int xsize, const unsigned int ysize,
			const std::vector<integral_t> wake);

public:
	/**
	 * @brief overloads HeritageMap::apply() to have a variable HeritageMap
	 */
	void apply();

private:
	/**
	 * @brief _wake (normalized) wake
	 */
	const std::vector<vfps::integral_t> _wake;

	/**
	 * @brief _wakeforce
	 */
	std::vector<vfps::meshaxis_t> _wakeforce;
};

} // namespace VFPS

#endif // KICKMAP_HPP
