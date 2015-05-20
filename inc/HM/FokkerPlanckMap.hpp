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

#ifndef FOKKERPLANCKMAP_HPP
#define FOKKERPLANCKMAP_HPP

#include "HeritageMap.hpp"

namespace vfps
{

/**
 * @brief The FokkerPlanckMap class
 *
 * Heritage information for the FokkerPlanckMap is one dimensional
 *
 * @todo Padding for faster memory access
 */
class FokkerPlanckMap : public HeritageMap
{
public:
	/**
	 * @brief The FPType enum holds ways to solve Fokker Planck equation
	 */
	enum class FPType : unsigned int {
		none=0,
		damping_only=1,
		diffusion_only=2,
		full=3
	};

public:
	/**
	 * @brief FokkerPlanckMap
	 * @param in
	 * @param out
	 * @param xsize
	 * @param ysize
	 * @param fpt
	 * @param e1 Marit: (deltat*2./(omegas*td))
	 */
	FokkerPlanckMap(PhaseSpace* in, PhaseSpace* out,
					const meshindex_t xsize, const meshindex_t ysize,
					FPType fpt, double e1);

	/**
	 * @brief apply custom apply method needed to handle one dimensional HM
	 */
	void apply();

private:
	/**
	 * @brief _meshxsize horizontal size of the meshaxis_t
	 *
	 * As this HeritageMap itself is one dimensional,
	 * it has to save the size of the actual mesh seperatly.
	 */
	const meshindex_t _meshxsize;
};

}

#endif // FOKKERPLANCKMAP_HPP
