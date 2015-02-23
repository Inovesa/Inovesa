/******************************************************************************/
/* Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   */
/* Copyright (c) 2007-2009: Peter Schregle (Fixed Point Math Library)         */
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

#ifndef ROTATIONMAP_HPP
#define ROTATIONMAP_HPP

#include "HeritageMap.hpp"

namespace vfps
{

class RotationMap : public HeritageMap
{
public:
	RotationMap(PhaseSpace* in, PhaseSpace* out,
				const size_t xsize, const size_t ysize,
				const meshaxis_t angle);
};

}

#endif // ROTATIONMAP_HPP
