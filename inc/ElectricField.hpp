#ifndef ELECTRICFIELD_HPP
#define ELECTRICFIELD_HPP

#include <algorithm>
#include <fftw3.h>
#include <sstream>

#include "defines.hpp"
#include "PhaseSpace.hpp"
#include "Ruler.hpp"
#include "Impedance.hpp"

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
	ElectricField(const PhaseSpace* phasespace, const Impedance* impedance);

	~ElectricField();

	inline csrpower_t* getData() const
		{ return _csrspectrum; }

	csrpower_t* updateCSRSpectrum();

private:
	fftwf_plan prepareFFT(unsigned int n, csrpower_t* in, impedance_t* out);

private:
	size_t _nmax;

	const integral_t* _bunchprofile;

	const size_t _bpmeshcells;

	csrpower_t* _csrspectrum;

	const Impedance* _impedance;

	integral_t* _bp_padded;

	integral_t* _bp_padded_fftw;

	impedance_t* _bp_fourier;

	fftwf_complex* _bp_fourier_fftw;

	fftwf_plan _ft_bunchprofile;

	const Ruler<meshaxis_t>* _spaceinfo;
};

} // namespace vfps

#endif // ELECTRICFIELD_HPP
