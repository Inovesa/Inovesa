/******************************************************************************/
/* Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   */
/* Copyright (c) 2014-2015: Patrik Schönfeldt                                 */
/*                                                                            */
/* This file is part of Inovesa.                                              */
/* Inovesa is free software: you can redistribute it and/or modify            */
/* it under the terms of the GNU General Public License as published by       */
/* the Free Software Foundation, either version 3 of the License, or          */
/* (at your option) any later version.                                        */
/*                                                                            */
/* Inovesa is distributed in the hope that it will be useful,                 */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           */
/******************************************************************************/

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
	ElectricField(PhaseSpace* phasespace, const Impedance* impedance);

	~ElectricField();

	inline csrpower_t* getData() const
		{ return _csrspectrum; }

	csrpower_t* updateCSRSpectrum();

private:
	fftwf_plan prepareFFT(unsigned int n, csrpower_t* in, impedance_t* out);

private:
	size_t _nmax;

	PhaseSpace* _phasespace;

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
