/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Algorithms   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
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
    Ruler(meshindex_t steps, meshaxis_t min, meshaxis_t max, double scale=0) :
        _steps(steps),
        _max(max),
        _min(min),
        _delta((max-min)/meshaxis_t(steps-1)),
        _scale(scale),
        _zerobin(((min+max)/(min-max)+1)*(steps-1)/2)
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
        _delta(other._delta),
        _scale(other._scale),
        _zerobin(other._zerobin)
    {
        meshaxis_t* meshaxis_tmp = new meshaxis_t[_steps];
        std::copy_n(other._data,_steps,meshaxis_tmp);
        _data = meshaxis_tmp;
    }

    ~Ruler()
    {
        delete [] _data;
    }

    inline const meshaxis_t* data() const
        { return _data; }

    inline const meshaxis_t max() const
        {return _max;}

    inline const meshaxis_t min() const
        {return _min;}

    inline double scale() const
        { return _scale; }

    inline meshindex_t steps() const
        {return _steps;}

    inline const meshaxis_t delta() const
        {return _delta;}

    inline const meshaxis_t size() const
        { return _max - _min; }

    inline const meshaxis_t zerobin() const
        { return _zerobin; }

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

    const double _scale;

    const meshaxis_t _zerobin;
};

}

#endif // RULER_HPP
