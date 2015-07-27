#include "Impedance.hpp"

vfps::Impedance::Impedance(vfps::Impedance::ImpedanceModel model, size_t n)
{
	_data.reserve(n);

	// according to Eq. 6.18 in Part. Acc. Vol 57, p 35 (Murpy et al.)
	constexpr impedance_t freespacecoeff = impedance_t(250.1,176.9);

	switch (model) {
	case ImpedanceModel::FreeSpace:
		for (size_t i=0; i<n; i++) {
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
}
