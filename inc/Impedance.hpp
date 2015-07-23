#ifndef IMPEDANCE_HPP
#define IMPEDANCE_HPP

#include <iostream>
#include <string>

#include "defines.hpp"

namespace vfps
{

class Impedance
{
public:
	Impedance(std::string datafile);

	inline const impedance_t* data() const
		{ return _data; }

	inline const impedance_t operator[](size_t n) const
		{ return _data[n]; }

private:
	const impedance_t* _data;
};

} // namespace vfps

#endif // IMPEDANCE_HPP
