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

namespace vfps
{

template <class data_t>
class Ruler
{
public:
	Ruler(size_t steps, data_t min, data_t max) :
        _steps(steps),
        _max(max),
        _min(min),
        _delta((max-min)/(steps-1))
    {
        if (min >= max) {
            throw std::invalid_argument("Tried to set up Ruler with min >= max.");
        }
        data_t* data_tmp = new data_t[_steps];

		for (size_t i=0; i<_steps; i++){
            data_tmp[i] = _min+(i*_delta);
        }

        _data = data_tmp;
    }

    Ruler(const Ruler& other) :
        _steps(other._steps),
        _max(other._max),
        _min(other._min),
        _delta(other._delta)
    {
        data_t* data_tmp = new data_t[_steps];
        std::copy_n(other._data,_steps,data_tmp);
        _data = data_tmp;
    }

    ~Ruler()
    {
        delete [] _data;
    }

    inline const data_t getMax() const
        {return _max;}

	inline const data_t getMin() const
        {return _min;}

	inline size_t getNSteps() const
        {return _steps;}

	inline const data_t getDelta() const
        {return _delta;}

	inline const data_t size() const
		{ return _max - _min; }

	inline const data_t& operator[](size_t d) const
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
    const data_t* _data;

	const size_t _steps;

    const data_t _max;

    const data_t _min;

	const data_t _delta;
};

}

#endif // RULER_HPP
