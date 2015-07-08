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

#ifndef RULER_HPP
#define RULER_HPP

#include <algorithm>
#include <stdexcept>

#include "defines.hpp"

namespace vfps
{

template <class meshaxis_t>
class Ruler
{
public:
	Ruler(meshindex_t steps, meshaxis_t min, meshaxis_t max) :
        _steps(steps),
        _max(max),
        _min(min),
		_delta((max-min)/meshaxis_t(steps-1))
    {
        if (min >= max) {
            throw std::invalid_argument("Tried to set up Ruler with min >= max.");
        }
		meshaxis_t* meshaxis_tmp = new meshaxis_t[_steps];

		for (meshindex_t i=0; i<_steps; i++){
			meshaxis_tmp[i] = _min+(meshaxis_t(i)*_delta);
        }

		_data = meshaxis_tmp;
    }

    Ruler(const Ruler& other) :
        _steps(other._steps),
        _max(other._max),
        _min(other._min),
        _delta(other._delta)
    {
		meshaxis_t* meshaxis_tmp = new meshaxis_t[_steps];
		std::copy_n(other._data,_steps,meshaxis_tmp);
		_data = meshaxis_tmp;
    }

    ~Ruler()
    {
        delete [] _data;
    }

	inline const meshaxis_t getMax() const
        {return _max;}

	inline const meshaxis_t getMin() const
        {return _min;}

	inline meshindex_t getNSteps() const
        {return _steps;}

	inline const meshaxis_t getDelta() const
        {return _delta;}

	inline const meshaxis_t size() const
		{ return _max - _min; }

	inline const meshaxis_t& operator[](meshindex_t d) const
        {return _data[d];}

    /**
     * @brief operator == compares grids
     * @param other grid to compare
     * @return true (same dimensions) or false (different dimensions)
     */
	bool operator==(const Ruler& other) const
    {
        if (_min == other._min && _max == other._max && _steps == other._steps){
            return true;
        } else {
            return false;
        }
    }

protected:
	const meshaxis_t* _data;

	const meshindex_t _steps;

	const meshaxis_t _max;

	const meshaxis_t _min;

	const meshaxis_t _delta;
};

}

#endif // RULER_HPP
