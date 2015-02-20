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
    Ruler(unsigned int steps, data_t min, data_t max) :
        _steps(steps),
        _max(max),
        _min(min),
        _delta((max-min)/(steps-1))
    {
        if (min >= max) {
            throw std::invalid_argument("Tried to set up Ruler with min >= max.");
        }
        data_t* data_tmp = new data_t[_steps];

        for (unsigned int i=0; i<_steps; i++){
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

	inline unsigned int getNSteps() const
        {return _steps;}

    inline const data_t getDelta() const
        {return _delta;}

	inline const data_t size() const
		{ return _max - _min; }

    inline const data_t& operator[](unsigned int d) const
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

    const unsigned int _steps;

    const data_t _max;

    const data_t _min;

    const data_t _delta;
};

}

#endif // RULER_HPP
