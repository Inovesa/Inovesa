#include "ElectricField.hpp"

vfps::ElectricField::ElectricField(const PhaseSpace* phasespace) :
	_nmax(phasespace->nMeshCells(0)),
	_bunchprofile(phasespace->getProjection(0)),
	_phasespace(phasespace),
	_spaceinfo(phasespace->getRuler(0))
{
	_bp_padded_fftw = fftw_alloc_real(2*_nmax);
	_bp_padded = reinterpret_cast<meshdata_t*>(_bp_padded_fftw);
	std::fill_n(&_bp_padded[_nmax],_nmax,integral_t(0)); //zero-padding

	_bp_fourier_fftw = fftw_alloc_complex(2*_nmax);
	_bp_fourier = reinterpret_cast<impedance_t*>(_bp_fourier_fftw);

	_ft_bunchprofile = prepareFFT(2*_nmax,_bp_padded,_bp_fourier);
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
	//FFT charge density
	fftw_execute(_ft_bunchprofile);

	for (unsigned int i=0; i<_nmax; i++) {
		// norm = squared magnitude
		_csrspectrum[i] = (_impedance[i]*std::norm(_bp_fourier[i])).real();
	}

	return _csrspectrum;
}

fftw_plan vfps::ElectricField::prepareFFT(	unsigned int n, csrpower_t* in,
											impedance_t* out)
{
	fftw_plan plan = nullptr;

	std::stringstream wisdomfile;
	wisdomfile << "wisdom_r2c_" << n << ".fftw";
	// use wisdomfile, if it exists
	if (fftw_import_wisdom_from_filename(wisdomfile.str().c_str()) != 0) {
		plan = fftw_plan_dft_r2c_1d(n,in,reinterpret_cast<fftw_complex*>(out),
									FFTW_WISDOM_ONLY|FFTW_PATIENT);
	}
	// if there was no wisdom (no or bad file), create some
	if (plan == nullptr) {
		plan = fftw_plan_dft_r2c_1d(n,in,reinterpret_cast<fftw_complex*>(out),
									FFTW_PATIENT);
		fftw_export_wisdom_to_filename(wisdomfile.str().c_str());
		Display::printText("Created some wisdom at "+wisdomfile.str());
	}
	return plan;
}
