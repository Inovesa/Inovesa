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

#ifndef WAKEKICKMAP_HPP
#define WAKEKICKMAP_HPP

#include <array>

#include "defines.hpp"
#include "KickMap.hpp"

using std::modf;

namespace vfps
{

/**
 * @brief The KickMap class offers an option for one-dimensional kicks
 *
 * @todo implement option to have 1D interpolation with HeritageMap
 * @todo change to use 1D HeritageMap
 */
class WakeKickMap : public KickMap
{
public:
	/**
	 * @brief WakeKickMap
	 * @param in
	 * @param out
	 * @param xsize
	 * @param ysize
	 * @param wakefunction from -xsize to xsize-1, normalized in a way
	 *	      that plain multiplication with density gives the force
	 */
	WakeKickMap(PhaseSpace* in, PhaseSpace* out,
				const meshindex_t xsize, const meshindex_t ysize,
				const std::pair<meshindex_t,double>* wake,
				const size_t wakesize, const InterpolationType it);

	~WakeKickMap();

public:
	/**
	 * @brief overloads HeritageMap::apply() to have a variable HeritageMap
	 */
	void apply();

private:
	/**
	 * @brief _wakefunktion (normalized single particle) wake,
	 *		  sampled at 2*xsize positions [-xsize:+xsize]
	 */
	meshaxis_t* const _wakefunction;
};

} // namespace VFPS

#endif // WAKEKICKMAP_HPP
