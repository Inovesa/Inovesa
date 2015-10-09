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
	 * @brief ElectricField minimal constructor, will not offer wake function
	 * @param phasespace this electric field is assigned to
	 * @param impedance to use for electric field calculation
	 */
	ElectricField(PhaseSpace* phasespace, const Impedance* impedance);

	/**
	 * @brief ElectricField
	 * @param phasespace this electric field is assigned to
	 * @param impedance to use for electric field calculation
	 * @param Ib bunch current [A]
	 * @param bl rms bunch length [m]
	 * @param E0 beam energy [eV]
	 * @param sigmaE normalized energy spread [1]
	 * @param fs synchrotron frequency [Hz]
	 * @param frev revolution frequency [Hz]
	 * @param dt time step [s]
	 * @param rbend bending radius [m]
	 */
	ElectricField(PhaseSpace* phasespace, const Impedance* impedance,
				  const double Ib, const double bl,
				  const double E0, const double sigmaE,
				  const double fs, const double frev,
				  const double dt, const double rbend);

	~ElectricField();

	inline csrpower_t* getData() const
		{ return _csrspectrum; }

	inline const Impedance* getImpedance() const
		{ return _impedance; }

	inline size_t getNMax() const
		{ return _nmax; }

	inline const Ruler<meshaxis_t>* getRuler() const
		{ return &_axis; }

	csrpower_t* updateCSRSpectrum();

	meshaxis_t* getWakefunction() const
		{ return _wakefunction; }

private:
	enum class fft_direction : uint_fast8_t {
		forward, backward
	};

	fftwf_plan prepareFFT(size_t n, csrpower_t* in, impedance_t* out);

	fftwf_plan prepareFFT(size_t n, impedance_t* in, impedance_t* out,
						  fft_direction direction);

private:
	const Ruler<meshaxis_t> _axis;

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

	meshaxis_t* _wakefunction;
};

} // namespace vfps

#endif // ELECTRICFIELD_HPP
