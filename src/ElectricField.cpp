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

#include "ElectricField.hpp"

vfps::ElectricField::ElectricField(PhaseSpace* phasespace,
								   const Impedance* impedance) :
	_nmax(impedance->maxN()),
	_phasespace(phasespace),
	_bpmeshcells(phasespace->nMeshCells(0)),
	_csrspectrum(new csrpower_t[_nmax]),
	_impedance(impedance),
	_spaceinfo(phasespace->getRuler(0))
{
	_bp_padded_fftw = fftwf_alloc_real(2*_nmax);
	_bp_padded = reinterpret_cast<meshdata_t*>(_bp_padded_fftw);

	//zero-padding
	std::fill_n(&_bp_padded[_bpmeshcells],2*_nmax-_bpmeshcells,integral_t(0));

	_bp_fourier_fftw = fftwf_alloc_complex(2*_nmax);
	_bp_fourier = reinterpret_cast<impedance_t*>(_bp_fourier_fftw);

	_ft_bunchprofile = prepareFFT(2*_nmax,_bp_padded,_bp_fourier);
}

vfps::ElectricField::~ElectricField()
{
	delete [] _csrspectrum;

	fftwf_free(_bp_padded_fftw);
	fftwf_free(_bp_fourier_fftw);
	fftwf_destroy_plan(_ft_bunchprofile);
	fftwf_cleanup();
}

vfps::csrpower_t* vfps::ElectricField::updateCSRSpectrum()
{
	std::copy_n(_phasespace->projectionToX(),
				_spaceinfo->getNSteps(),
				_bp_padded);

	//FFT charge density
	fftwf_execute(_ft_bunchprofile);

	for (unsigned int i=0; i<_nmax; i++) {
		// norm = squared magnitude
		_csrspectrum[i] = ((*_impedance)[i]*std::norm(_bp_fourier[i])).real();
	}

	return _csrspectrum;
}

fftwf_plan vfps::ElectricField::prepareFFT(	unsigned int n, csrpower_t* in,
											impedance_t* out)
{
	fftwf_plan plan = nullptr;

	std::stringstream wisdomfile;
	// get ready to save BunchCharge
	if (std::is_same<vfps::csrpower_t,float>::value) {
		wisdomfile << "wisdom_r2c32_" << n << ".fftw";
	} else if (std::is_same<vfps::csrpower_t,double>::value) {
		wisdomfile << "wisdom_r2c64_" << n << ".fftw";
	}
	// use wisdomfile, if it exists
	if (fftwf_import_wisdom_from_filename(wisdomfile.str().c_str()) != 0) {
		plan = fftwf_plan_dft_r2c_1d(n,in,reinterpret_cast<fftwf_complex*>(out),
									FFTW_WISDOM_ONLY|FFTW_PATIENT);
	}
	// if there was no wisdom (no or bad file), create some
	if (plan == nullptr) {
		plan = fftwf_plan_dft_r2c_1d(n,in,reinterpret_cast<fftwf_complex*>(out),
									FFTW_PATIENT);
		fftwf_export_wisdom_to_filename(wisdomfile.str().c_str());
		Display::printText("Created some wisdom at "+wisdomfile.str());
	}
	return plan;
}
