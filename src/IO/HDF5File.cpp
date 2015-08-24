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

#include "IO/HDF5File.hpp"

vfps::HDF5File::HDF5File(std::string fname,
						 const PhaseSpace* ps,
						 const size_t maxn) :
	file( nullptr ),
	fname( fname ),
	bc_dims( 0 ),
	bp_dims( {{ 0, ps->nMeshCells(0) }} ),
	csr_dims( {{ 0, maxn }} ),
	maxn(maxn),
	ps_dims( {{ 0, ps->nMeshCells(0), ps->nMeshCells(1) }} ),
	ps_size( ps->nMeshCells(0) )
{
	file = new H5::H5File(fname,H5F_ACC_TRUNC);

	file->createGroup("Info");

	// save Values of Axis
	if (std::is_same<vfps::meshaxis_t,float>::value) {
		psa_datatype = H5::PredType::IEEE_F32LE;
	} else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
		psa_datatype = H5::PredType::STD_I64LE;
	} else if (std::is_same<vfps::meshaxis_t,double>::value) {
		psa_datatype = H5::PredType::IEEE_F64LE;
	}

	const std::array<hsize_t,psa_rank> psa_dims
			= {{ 1, ps_size }};

	const std::array<hsize_t,psa_rank> psa_maxdims
			= {{ 1, ps_size }};

	psa0_dataspace = new H5::DataSpace(psa_rank,psa_dims.data(),
									   psa_maxdims.data());
	psa1_dataspace = new H5::DataSpace(psa_rank,psa_dims.data(),
									   psa_maxdims.data());


	const std::array<hsize_t,bp_rank> psa_chunkdims
			= {{1U,ps_size/8U}};
	psa_prop.setChunk(psa_rank,psa_chunkdims.data());
	psa_prop.setShuffle();
	psa_prop.setDeflate(compression);

	psa0_dataset = file->createDataSet("/Info/AxisValues_z",psa_datatype,
											*psa0_dataspace,psa_prop);

	psa0_dataset.write(ps->getRuler(0)->data(),psa_datatype);

	psa1_dataset = file->createDataSet("/Info/AxisValues_E",psa_datatype,
											*psa1_dataspace,psa_prop);

	psa1_dataset.write(ps->getRuler(0)->data(),psa_datatype);

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

	const hsize_t bc_chunkdims = 8;
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
			= {{1U,ps_size/8U}};
	bp_prop.setChunk(bp_rank,bp_chunkdims.data());
	bp_prop.setShuffle();
	bp_prop.setDeflate(compression);

	bp_dataset = file->createDataSet("/BunchProfile/data",bp_datatype,
											*bp_dataspace,bp_prop);

	// get ready to save CSR Spectrum
	file->createGroup("CSR-Spectrum");

	if (std::is_same<vfps::csrpower_t,float>::value) {
		csr_datatype = H5::PredType::IEEE_F32LE;
	} else if (std::is_same<vfps::csrpower_t,double>::value) {
		csr_datatype = H5::PredType::IEEE_F64LE;
	}

	const std::array<hsize_t,csr_rank> csr_maxdims
			= {{H5S_UNLIMITED,maxn}};

	csr_dataspace = new H5::DataSpace(csr_rank,csr_dims.data(),
									 csr_maxdims.data());


	const std::array<hsize_t,csr_rank> csr_chunkdims
			= {{1U,ps_size/8U}};
	csr_prop.setChunk(csr_rank,csr_chunkdims.data());
	csr_prop.setShuffle();
	csr_prop.setDeflate(compression);

	csr_dataset = file->createDataSet("/CSR-Spectrum/data",csr_datatype,
											*csr_dataspace,csr_prop);

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
			= {{1U,ps_size/8U,ps_size/8U}};
	ps_prop.setChunk(ps_rank,ps_chunkdims.data());
	ps_prop.setShuffle();
	ps_prop.setDeflate(compression);

	ps_dataset = file->createDataSet("/PhaseSpace/data",ps_datatype,
											*ps_dataspace,ps_prop);

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
	delete psa0_dataspace;
	delete psa1_dataspace;
	delete bc_dataspace;
	delete bp_dataspace;
	delete ps_dataspace;
}

void vfps::HDF5File::append(const ElectricField* ef)
{
	// append CSR Field
	std::array<hsize_t,csr_rank> csr_offset
			= {{csr_dims[0],0}};
	const std::array<hsize_t,csr_rank> csr_ext
			= {{1,maxn}};
	csr_dims[0]++;
	csr_dataset.extend(csr_dims.data());
	H5::DataSpace* filespace = new H5::DataSpace(csr_dataset.getSpace());
	filespace->selectHyperslab(H5S_SELECT_SET, csr_ext.data(), csr_offset.data());
	H5::DataSpace* memspace = new H5::DataSpace(csr_rank,csr_ext.data(),nullptr);
	csr_dataset.write(ef->getData(), csr_datatype,*memspace, *filespace);
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

	// append BunchCharge
	hsize_t bc_offset = bc_dims;
	const hsize_t bc_ext = 1;
	bc_dims++;
	bc_dataset.extend(&bc_dims);
	filespace = new H5::DataSpace(bc_dataset.getSpace());
	filespace->selectHyperslab(H5S_SELECT_SET, &bc_ext, &bc_offset);
	memspace = new H5::DataSpace(bc_rank,&bc_ext,nullptr);
	integral_t bunchcharge = ps->getCharge();
	bc_dataset.write(&bunchcharge, bc_datatype,*memspace, *filespace);
	delete memspace;
	delete filespace;
}
