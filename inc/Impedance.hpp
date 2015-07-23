#ifndef IMPEDANCE_HPP
#define IMPEDANCE_HPP

#include <fstream>
#include <string>
#include <vector>

#include "defines.hpp"

namespace vfps
{

class Impedance
{
public:
	Impedance(std::string datafile);

	inline const impedance_t* data() const
		{ return _data.data(); }

	inline const impedance_t operator[](size_t n) const
		{ return _data[n]; }

	inline size_t maxN() const
		{ return _data.size(); }

private:
	std::vector<impedance_t> _data;
};

} // namespace vfps

#endif // IMPEDANCE_HPP
