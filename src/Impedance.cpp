#include "Impedance.hpp"

vfps::Impedance::Impedance(std::string datafile)
{
	std::ifstream is(datafile);
	size_t lineno;
	double real;
	double imag;

	while(is.good()) {
		is >> lineno >> real >> imag;
		_data.push_back(impedance_t(real,imag));
	}
}
