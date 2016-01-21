/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   *
 * Copyright (c) 2014-2015: Patrik Schönfeldt                                 *
 *                                                                            *
 * This file is part of Inovesa.                                              *
 * Inovesa is free software: you can redistribute it and/or modify            *
 * it under the terms of the GNU General Public License as published by       *
 * the Free Software Foundation, either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * Inovesa is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU General Public License for more details.                               *
 *                                                                            *
 * You should have received a copy of the GNU General Public License          *
 * along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           *
 ******************************************************************************/

#ifdef INOVESA_USE_HDF5
#include "IO/HDF5File.hpp"

vfps::HDF5File::HDF5File(const std::string fname,
                         const PhaseSpace* ps,
                         const ElectricField* ef,
                         const Impedance* imp,
                         const WakeFunctionMap* wfm) :
    file( nullptr ),
    fname( fname ),
    ta_dims( 0 ),
    bc_dims( 0 ),
    bp_dims( {{ 0, ps->nMeshCells(0) }} ),
    bl_dims( 0 ),
    wp_dims( {{ 0, ps->nMeshCells(0) }} ),
    csr_dims( {{ 0, ef->getNMax() }} ),
    maxn( ef->getNMax() ),
    csrp_dims( 0 ),
    ps_dims( {{ 0, ps->nMeshCells(0), ps->nMeshCells(1) }} ),
    ps_size( ps->nMeshCells(0) ),
    imp_size( imp->maxN() ),
    wf_size( 2*ps_size )
{
    file = new H5::H5File(fname,H5F_ACC_TRUNC);

	file->createGroup("Info");
	// save Values of Phase Space Axis
	if (std::is_same<vfps::meshaxis_t,float>::value) {
		axps_datatype = H5::PredType::IEEE_F32LE;
	} else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
		axps_datatype = H5::PredType::STD_I64LE;
	} else if (std::is_same<vfps::meshaxis_t,double>::value) {
		axps_datatype = H5::PredType::IEEE_F64LE;
	}

	const std::array<hsize_t,axps_rank> psa_dims
			= {{ 1, ps_size }};
	const std::array<hsize_t,axps_rank> psa_maxdims
			= {{ 1, ps_size }};

	ax0ps_dataspace = new H5::DataSpace(axps_rank,psa_dims.data(),
									   psa_maxdims.data());
	ax1ps_dataspace = new H5::DataSpace(axps_rank,psa_dims.data(),
									   psa_maxdims.data());


	const std::array<hsize_t,axps_rank> psa_chunkdims
			= {{1U,std::min(2048U,ps_size)}};
	axps_prop.setChunk(axps_rank,psa_chunkdims.data());
	axps_prop.setShuffle();
	axps_prop.setDeflate(compression);

	ax0ps_dataset = file->createDataSet("/Info/AxisValues_z",axps_datatype,
											*ax0ps_dataspace,axps_prop);
	const double ax0scale = ps->getScale(0);
	ax0ps_dataset.createAttribute("Scale",H5::PredType::IEEE_F64LE,
		H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax0scale);
	ax1ps_dataset = file->createDataSet("/Info/AxisValues_E",axps_datatype,
											*ax1ps_dataspace,axps_prop);
	const double ax1scale = ps->getScale(1);
	ax1ps_dataset.createAttribute("Scale",H5::PredType::IEEE_F64LE,
		H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax1scale);
	ax0ps_dataset.write(ps->getRuler(0)->data(),axps_datatype);
	ax1ps_dataset.write(ps->getRuler(0)->data(),axps_datatype);

	// save Values of Frequency Axis
	if (std::is_same<vfps::meshaxis_t,float>::value) {
		axfreq_datatype = H5::PredType::IEEE_F32LE;
	} else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
		axfreq_datatype = H5::PredType::STD_I64LE;
	} else if (std::is_same<vfps::meshaxis_t,double>::value) {
		axfreq_datatype = H5::PredType::IEEE_F64LE;
	}

	const std::array<hsize_t,axfreq_rank> axfreq_dims
			= {{ 1, maxn }};
	const std::array<hsize_t,axfreq_rank> axfreq_maxdims
			= {{ 1, maxn }};

    axfreq_dataspace = new H5::DataSpace(axfreq_rank,axfreq_dims.data(),
                                         axfreq_maxdims.data());


	const std::array<hsize_t,bp_rank> axfreq_chunkdims
			= {{1U,std::min(2048U,ps_size)}};
	axfreq_prop.setChunk(axps_rank,axfreq_chunkdims.data());
	axfreq_prop.setShuffle();
	axfreq_prop.setDeflate(compression);

    axfreq_dataset = file->createDataSet("/Info/AxisValues_f",axfreq_datatype,
                                         *axfreq_dataspace,axfreq_prop);
    const double axfreqscale = imp->getRuler()->scale();
    axfreq_dataset.createAttribute("Scale",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&axfreqscale);
    axfreq_dataset.write(imp->getRuler()->data(),axfreq_datatype);


    // get ready to save TimeAxis
    ta_datatype = H5::PredType::IEEE_F64LE;

    hsize_t ta_maxdims = H5S_UNLIMITED;

    ta_dataspace = new H5::DataSpace(ta_rank,&ta_dims,&ta_maxdims);


    const hsize_t ta_chunkdims = 256;
    ta_prop.setChunk(ta_rank,&ta_chunkdims);
    ta_prop.setShuffle();
    ta_prop.setDeflate(compression);

    ta_dataset = file->createDataSet("/Info/AxisValues_t",ta_datatype,
                                     *ta_dataspace,ta_prop);

	// get ready to save BunchCharge
	file->createGroup("BunchCharge");
	if (std::is_same<vfps::integral_t,float>::value) {
		bc_datatype = H5::PredType::IEEE_F32LE;
	} else if (std::is_same<vfps::integral_t,fixp64>::value) {
		bc_datatype = H5::PredType::STD_I64LE;
	} else if (std::is_same<vfps::integral_t,double>::value) {
		bc_datatype = H5::PredType::IEEE_F64LE;
	}

	hsize_t bc_maxdims = H5S_UNLIMITED;

	bc_dataspace = new H5::DataSpace(bc_rank,&bc_dims,&bc_maxdims);

	const hsize_t bc_chunkdims = 256;
	bc_prop.setChunk(bc_rank,&bc_chunkdims);
	bc_prop.setShuffle();
	bc_prop.setDeflate(compression);

	bc_dataset = file->createDataSet("/BunchCharge/data",bc_datatype,
											*bc_dataspace,bc_prop);

	// get ready to save BunchProfiles
	file->createGroup("BunchProfile");
	file->link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/BunchProfile/axis0" );

	if (std::is_same<vfps::integral_t,float>::value) {
		bp_datatype = H5::PredType::IEEE_F32LE;
	} else if (std::is_same<vfps::integral_t,fixp64>::value) {
		bp_datatype = H5::PredType::STD_I64LE;
	} else if (std::is_same<vfps::integral_t,double>::value) {
		bp_datatype = H5::PredType::IEEE_F64LE;
	}

	const std::array<hsize_t,bp_rank> bp_maxdims
			= {{H5S_UNLIMITED,ps_size}};

	bp_dataspace = new H5::DataSpace(bp_rank,bp_dims.data(),bp_maxdims.data());

    const std::array<hsize_t,bp_rank> bp_chunkdims
        = {{64U,std::min(2048U,ps_size)}};
    bp_prop.setChunk(bp_rank,bp_chunkdims.data());
    bp_prop.setShuffle();
    bp_prop.setDeflate(compression);

	bp_dataset = file->createDataSet("/BunchProfile/data",bp_datatype,
					 *bp_dataspace,bp_prop);

    // get ready to save BunchLength
    file->createGroup("BunchLength");
    if (std::is_same<vfps::meshaxis_t,float>::value) {
            bl_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            bl_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            bl_datatype = H5::PredType::IEEE_F64LE;
    }

    hsize_t bl_maxdims = H5S_UNLIMITED;

    bl_dataspace = new H5::DataSpace(bl_rank,&bl_dims,&bl_maxdims);

    const hsize_t bl_chunkdims = 256;
    bl_prop.setChunk(bl_rank,&bl_chunkdims);
    bl_prop.setShuffle();
    bl_prop.setDeflate(compression);

    bl_dataset = file->createDataSet("/BunchLength/data",bl_datatype,
                                     *bl_dataspace,bl_prop);

    // get ready to save WakePotential
    file->createGroup("WakePotential");
    file->link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/WakePotential/axis0" );

    if (std::is_same<vfps::meshaxis_t,float>::value) {
        wp_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
        wp_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
        wp_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,wp_rank> wp_maxdims
            = {{H5S_UNLIMITED,ps_size}};

    wp_dataspace = new H5::DataSpace(wp_rank,wp_dims.data(),wp_maxdims.data());

    const std::array<hsize_t,wp_rank> wp_chunkdims
            = {{64U,std::min(2048U,ps_size)}};
    wp_prop.setChunk(wp_rank,wp_chunkdims.data());
    wp_prop.setShuffle();
    wp_prop.setDeflate(compression);

    wp_dataset = file->createDataSet("/WakePotential/data",wp_datatype,
                                     *wp_dataspace,wp_prop);

    // get ready to save CSR Spectrum
    file->createGroup("CSRSpectrum");
    file->link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/CSRSpectrum/axis0" );
    file->link(H5L_TYPE_SOFT, "/Info/AxisValues_f", "/CSRSpectrum/axis1" );

    if (std::is_same<vfps::csrpower_t,float>::value) {
        csr_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::csrpower_t,double>::value) {
        csr_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,csr_rank> csr_maxdims
        = {{H5S_UNLIMITED,csr_dims[1]}};

    csr_dataspace = new H5::DataSpace(csr_rank,csr_dims.data(),
                                      csr_maxdims.data());


    const std::array<hsize_t,csr_rank> csr_chunkdims
        = {{64U,std::min(2048U,ps_size)}};
    csr_prop.setChunk(csr_rank,csr_chunkdims.data());
    csr_prop.setShuffle();
    csr_prop.setDeflate(compression);

    csr_dataset = file->createDataSet("/CSRSpectrum/data",csr_datatype,
                                      *csr_dataspace,csr_prop);

    // get ready to save CSR Power
    file->createGroup("CSRPower");
    file->link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/CSRPower/axis0" );

    if (std::is_same<vfps::csrpower_t,float>::value) {
        csrp_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::csrpower_t,double>::value) {
        csrp_datatype = H5::PredType::IEEE_F64LE;
    }

    const hsize_t csrp_maxdims = H5S_UNLIMITED;

    csrp_dataspace = new H5::DataSpace(csrp_rank,&csrp_dims,&csrp_maxdims);


    const hsize_t csrp_chunkdims = std::min(2048U,ps_size);
    csrp_prop.setChunk(csrp_rank,&csrp_chunkdims);
    csrp_prop.setShuffle();
    csrp_prop.setDeflate(compression);

    csrp_dataset = file->createDataSet("/CSRPower/data",csrp_datatype,
                                      *csrp_dataspace,csrp_prop);

	// get ready to save PhaseSpace
	file->createGroup("PhaseSpace");
	file->link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/PhaseSpace/axis0" );
	file->link(H5L_TYPE_SOFT, "/Info/AxisValues_E", "/PhaseSpace/axis1" );

	if (std::is_same<vfps::meshdata_t,float>::value) {
		ps_datatype = H5::PredType::IEEE_F32LE;
#if FXP_FRACPART < 31
	} else if (std::is_same<vfps::meshdata_t,fixp32>::value) {
		ps_datatype = H5::PredType::STD_I32LE;
#endif
	} else if (std::is_same<vfps::meshdata_t,fixp64>::value) {
		ps_datatype = H5::PredType::STD_I64LE;
	}else if (std::is_same<vfps::meshdata_t,double>::value) {
		ps_datatype = H5::PredType::IEEE_F64LE;
	}

	const std::array<hsize_t,ps_rank> ps_maxdims
			= {{H5S_UNLIMITED,ps_size,ps_size}};

	ps_dataspace = new H5::DataSpace(ps_rank,ps_dims.data(),ps_maxdims.data());


    const std::array<hsize_t,ps_rank> ps_chunkdims
        = {{64U,std::min(64U,ps_size),std::min(64U,ps_size)}};
    ps_prop.setChunk(ps_rank,ps_chunkdims.data());
    ps_prop.setShuffle();
    ps_prop.setDeflate(compression);

	ps_dataset = file->createDataSet("/PhaseSpace/data",ps_datatype,
											*ps_dataspace,ps_prop);

    // save Impedance
    file->createGroup("Impedance");
    file->link(H5L_TYPE_SOFT, "/Info/AxisValues_f", "/Impedance/axis0" );

    if (std::is_same<vfps::impedance_t,std::complex<float>>::value) {
        imp_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::impedance_t,std::complex<fixp64>>::value) {
        imp_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::impedance_t,std::complex<double>>::value) {
        imp_datatype = H5::PredType::IEEE_F64LE;
    }

    imp_dataspace = new H5::DataSpace(imp_rank,&imp_size,&imp_size);

    const hsize_t imp_chunkdims = std::min(hsize_t(4096),imp_size);
    imp_prop.setChunk(imp_rank,&imp_chunkdims);
    imp_prop.setShuffle();
    imp_prop.setDeflate(compression);

    file->createGroup("Impedance/data");
    imp_dataset_real = file->createDataSet("/Impedance/data/real",imp_datatype,
                                           *imp_dataspace,imp_prop);
    imp_dataset_imag = file->createDataSet("/Impedance/data/imag",imp_datatype,
                                           *imp_dataspace,imp_prop);

    std::vector<csrpower_t> imp_real;
    std::vector<csrpower_t> imp_imag;
    imp_real.reserve(imp_size);
    imp_imag.reserve(imp_size);
    for (const impedance_t& z : imp->impedance()) {
        imp_real.push_back(z.real());
        imp_imag.push_back(z.imag());
    }
    imp_dataset_real.write(imp_real.data(),imp_datatype);
    imp_dataset_imag.write(imp_imag.data(),imp_datatype);

    if (wfm != nullptr ) {
    // save Wake Function
    file->createGroup("WakeFunction");

    if (std::is_same<vfps::meshaxis_t,float>::value) {
        wf_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
        wf_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
        wf_datatype = H5::PredType::IEEE_F64LE;
    }

    wf_dataspace = new H5::DataSpace(wf_rank,&wf_size,&wf_size);

    const hsize_t wf_chunkdims = std::min(hsize_t(4096),wf_size);
    wf_prop.setChunk(wf_rank,&wf_chunkdims);
    wf_prop.setShuffle();
    wf_prop.setDeflate(compression);

    wf_dataset = file->createDataSet("/WakeFunction/data",wf_datatype,
                                     *wf_dataspace,wf_prop);

    wf_dataset.write(wfm->getWakeFunction(),wf_datatype);
    }

    // save Inovesa version
        std::array<hsize_t,1> version_dims {{3}};
        H5::DataSpace version_dspace(1,version_dims.data(),version_dims.data());
        H5::DataSet version_dset = file->createDataSet
                        ("/Info/INOVESA_v", H5::PredType::STD_I32LE,version_dspace);
        std::array<int32_t,3> version {{INOVESA_VERSION_RELEASE,
                                                                        INOVESA_VERSION_MINOR,
                                                                        INOVESA_VERSION_FIX}};
        version_dset.write(version.data(),H5::PredType::NATIVE_INT);
}

vfps::HDF5File::~HDF5File()
{
	delete file;
	delete axfreq_dataspace;
	delete ax0ps_dataspace;
	delete ax1ps_dataspace;
	delete bc_dataspace;
	delete bp_dataspace;
	delete ps_dataspace;
}

void vfps::HDF5File::append(const ElectricField* ef)
{
    H5::DataSpace* filespace;
    H5::DataSpace* memspace;

	// append CSR Field
	std::array<hsize_t,csr_rank> csr_offset
			= {{csr_dims[0],0}};
	const std::array<hsize_t,csr_rank> csr_ext
			= {{1,csr_dims[1]}};
	csr_dims[0]++;
	csr_dataset.extend(csr_dims.data());
	filespace = new H5::DataSpace(csr_dataset.getSpace());
	filespace->selectHyperslab(H5S_SELECT_SET, csr_ext.data(), csr_offset.data());
	memspace = new H5::DataSpace(csr_rank,csr_ext.data(),nullptr);
	csr_dataset.write(ef->getCSRSpectrum(), csr_datatype,*memspace, *filespace);
	delete memspace;
	delete filespace;

    // append CSR Power
    hsize_t csrp_offset = csrp_dims;
    const hsize_t csrp_ext = 1;
    csrp_dims++;
    csrp_dataset.extend(&csrp_dims);
    filespace = new H5::DataSpace(csrp_dataset.getSpace());
    filespace->selectHyperslab(H5S_SELECT_SET, &csrp_ext, &csrp_offset);
    memspace = new H5::DataSpace(csrp_rank,&csrp_ext,nullptr);
    csrpower_t csrpower = ef->getCSRPower();
    csrp_dataset.write(&csrpower, csrp_datatype,*memspace, *filespace);
    delete memspace;
    delete filespace;
}

void vfps::HDF5File::append(const PhaseSpace* ps)
{
	// append PhaseSpace
	std::array<hsize_t,ps_rank> ps_offset
			= {{ps_dims[0],0,0}};
	const std::array<hsize_t,ps_rank> ps_ext
			= {{1,ps_size,ps_size}};
	ps_dims[0]++;
	ps_dataset.extend(ps_dims.data());
	H5::DataSpace* filespace = new H5::DataSpace(ps_dataset.getSpace());
	H5::DataSpace* memspace = new H5::DataSpace(ps_rank,ps_ext.data(),nullptr);
	filespace->selectHyperslab(H5S_SELECT_SET, ps_ext.data(), ps_offset.data());
	ps_dataset.write(ps->getData(), ps_datatype, *memspace, *filespace);
	delete memspace;
	delete filespace;

	// append BunchProfile
	std::array<hsize_t,bp_rank> bp_offset
			= {{bp_dims[0],0}};
	const std::array<hsize_t,bp_rank> bp_ext
			= {{1,ps_size}};
	bp_dims[0]++;
	bp_dataset.extend(bp_dims.data());
	filespace = new H5::DataSpace(bp_dataset.getSpace());
	filespace->selectHyperslab(H5S_SELECT_SET, bp_ext.data(), bp_offset.data());
	memspace = new H5::DataSpace(bp_rank,bp_ext.data(),nullptr);
	bp_dataset.write(ps->getProjection(0), bp_datatype,*memspace, *filespace);
	delete memspace;
	delete filespace;

    // append Length
    hsize_t bl_offset = bl_dims;
    const hsize_t bl_ext = 1;
    bl_dims++;
    bl_dataset.extend(&bl_dims);
    filespace = new H5::DataSpace(bl_dataset.getSpace());
    filespace->selectHyperslab(H5S_SELECT_SET, &bl_ext, &bl_offset);
    memspace = new H5::DataSpace(bl_rank,&bl_ext,nullptr);
    meshaxis_t bunchlength = ps->getMoment(0,1);
    bl_dataset.write(&bunchlength, bl_datatype,*memspace, *filespace);
    delete memspace;
    delete filespace;

	// append BunchCharge
	hsize_t bc_offset = bc_dims;
	const hsize_t bc_ext = 1;
	bc_dims++;
	bc_dataset.extend(&bc_dims);
	filespace = new H5::DataSpace(bc_dataset.getSpace());
	filespace->selectHyperslab(H5S_SELECT_SET, &bc_ext, &bc_offset);
	memspace = new H5::DataSpace(bc_rank,&bc_ext,nullptr);
	integral_t bunchcharge = ps->getIntegral();
	bc_dataset.write(&bunchcharge, bc_datatype,*memspace, *filespace);
	delete memspace;
	delete filespace;
}

void vfps::HDF5File::append(const WakeKickMap* wkm)
{
    // append WakePotential
    std::array<hsize_t,wp_rank> wp_offset
            = {{wp_dims[0],0}};
    const std::array<hsize_t,wp_rank> wp_ext
            = {{1,ps_size}};
    wp_dims[0]++;
    wp_dataset.extend(wp_dims.data());
    H5::DataSpace* filespace = new H5::DataSpace(wp_dataset.getSpace());
    filespace->selectHyperslab(H5S_SELECT_SET, wp_ext.data(), wp_offset.data());
    H5::DataSpace* memspace = new H5::DataSpace(wp_rank,wp_ext.data(),nullptr);
    wp_dataset.write(wkm->getForce(), wp_datatype,*memspace, *filespace);
    delete memspace;
    delete filespace;
}

void vfps::HDF5File::timeStep(const double t)
{
    // append to TimeAxis
    hsize_t ta_offset = ta_dims;
    const hsize_t ta_ext = 1;
    ta_dims++;
    ta_dataset.extend(&ta_dims);
    H5::DataSpace* filespace;
    H5::DataSpace* memspace;
    filespace = new H5::DataSpace(ta_dataset.getSpace());
    filespace->selectHyperslab(H5S_SELECT_SET, &ta_ext, &ta_offset);
    memspace = new H5::DataSpace(ta_rank,&ta_ext,nullptr);
    ta_dataset.write(&t, ta_datatype,*memspace, *filespace);
    delete memspace;
    delete filespace;
}

#endif // INOVESA_USE_HDF5
