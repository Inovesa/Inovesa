#include "ElectricField.hpp"

vfps::ElectricField::ElectricField(const PhaseSpace* phasespace) :
	_nmax(phasespace->nMeshCells(0)),
	_bunchprofile(phasespace->getProjection(0)),
	_phasespace(phasespace),
	_spaceinfo(phasespace->getRuler(0))
{
	_bp_padded_fftw = fftw_alloc_real(2*_nmax);
	_bp_padded = reinterpret_cast<meshdata_t*>(_bp_padded_fftw);
	_bp_fourier_fftw = fftw_alloc_complex(2*_nmax);
	_bp_fourier = reinterpret_cast<impedance_t*>(_bp_padded_fftw);

	_ft_bunchprofile = fftw_plan_dft_r2c_1d(2*_nmax,_bp_padded_fftw,
											_bp_fourier_fftw,FFTW_ESTIMATE);
}

vfps::ElectricField::~ElectricField()
{
	fftw_free(_bp_padded_fftw);
	fftw_free(_bp_fourier_fftw);
	fftw_destroy_plan(_ft_bunchprofile);
	fftw_cleanup();
}

vfps::csrpower_t* vfps::ElectricField::updateCSRSpectrum()
{
	for (unsigned int i=0; i< _nmax; i++) {
		_bp_padded[i] = _bunchprofile[i];
	}
	for (unsigned int i= _nmax; i<2*_nmax; i++){ //zero-padding
		_bp_padded[i] = 0.0;
	}
	//FFT charge density
	fftw_execute(_ft_bunchprofile);

	for (unsigned int i=0; i<_nmax; i++) {
		// norm = squared magnitude
		_csrspectrum[i] = (_impedance[i]*std::norm(_bp_padded[i])).real();
	}

	return _csrspectrum;
}
