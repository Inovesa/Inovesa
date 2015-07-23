#ifndef ELECTRICFIELD_HPP
#define ELECTRICFIELD_HPP

#include <algorithm>
#include <fftw3.h>
#include <sstream>

#include "defines.hpp"
#include "PhaseSpace.hpp"
#include "Ruler.hpp"

#include "IO/Display.hpp"

namespace vfps {

class ElectricField
{
public:
	/**
	 * @brief ElectricField
	 *
	 * @param phasespace this electric field is assigned to
	 */
	ElectricField(const PhaseSpace* phasespace);

	~ElectricField();

	csrpower_t* updateCSRSpectrum();

private:
	fftw_plan prepareFFT(unsigned int n, csrpower_t* in, impedance_t* out);

private:
	size_t _nmax;

	const integral_t* _bunchprofile;

	csrpower_t* _csrspectrum;

	impedance_t* _impedance;

	integral_t* _bp_padded;

	integral_t* _bp_padded_fftw;

	impedance_t* _bp_fourier;

	fftw_complex* _bp_fourier_fftw;

	fftw_plan _ft_bunchprofile;

	const Ruler<meshaxis_t>* _spaceinfo;
};

} // namespace vfps

#endif // ELECTRICFIELD_HPP
