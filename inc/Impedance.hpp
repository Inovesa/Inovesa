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
	enum class ImpedanceModel : uint_fast8_t {
		FreeSpace = 0
	};

public:
	/**
	 * @brief Impedance
	 * @param model
	 * @param n
	 * @param roundup round to next 2^N
	 */
	Impedance(ImpedanceModel model, size_t n, bool roundup=true);

	Impedance(std::string datafile);

	inline const impedance_t* data() const
		{ return _data.data(); }

	inline const impedance_t operator[](size_t n) const
		{ return _data[n]; }

	inline size_t maxN() const
		{ return _data.size(); }

private:
	std::vector<impedance_t> _data;

	size_t _nmax;

	/**
	 * @brief upper_power_of_two
	 * @param v
	 * @return
	 *
	 * see http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
	 */
	uint64_t upper_power_of_two(uint64_t v);


};

} // namespace vfps

#endif // IMPEDANCE_HPP
