#include "Impedance.hpp"

vfps::Impedance::Impedance(vfps::Impedance::ImpedanceModel model,
						   size_t n, bool roundup)
{
	if (roundup) {
		   _nmax = upper_power_of_two(n);
	} else {
		_nmax = n;
	}
	_data.reserve(_nmax);

	// according to Eq. 6.18 in Part. Acc. Vol 57, p 35 (Murpy et al.)
	constexpr impedance_t freespacecoeff = impedance_t(250.1,176.9);

	switch (model) {
	case ImpedanceModel::FreeSpace:
		for (size_t i=0; i<_nmax; i++) {
			_data.push_back(freespacecoeff*std::pow(csrpower_t(n),
													csrpower_t(1.0/3.0)));
		}
	break;
	}
}

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
	_nmax = _data.size();
}

uint64_t vfps::Impedance::upper_power_of_two(uint64_t v)
{
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	v |= v >> 32;
	v++;
	return v;
}
