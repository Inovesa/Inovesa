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

vfps::HDF5File::HDF5File(std::string fname) :
	file( nullptr ),
	fname( fname ),
	bp_dims( {{ ps_xsize,1 }}),
	bp_name( "BunchProfile" ),
	ps_dims( {{ ps_xsize,ps_ysize,1}} ),
	ps_name( "PhaseSpace" )
{
	file = new H5::H5File(fname,H5F_ACC_TRUNC);

	// get ready to save BunchProfiles
	if (std::is_same<vfps::integral_t,float>::value) {
		bp_datatype = H5::PredType::IEEE_F32LE;
	} else {
		bp_datatype = H5::PredType::STD_I64LE;
	}

	static constexpr std::array<hsize_t,bp_rank> bp_maxdims
			= {{ps_xsize,H5S_UNLIMITED}};

	bp_dataspace = new H5::DataSpace(bp_rank,bp_dims.data(),bp_maxdims.data());


	static constexpr std::array<hsize_t,bp_rank> bp_chunkdims
			= {{ps_xsize/8,1}};
	bp_prop.setChunk(bp_rank,bp_chunkdims.data());
	bp_prop.setShuffle();
	bp_prop.setDeflate(compression);

	bp_dataset = new H5::DataSet(
						file->createDataSet(bp_name,bp_datatype,
											*bp_dataspace,bp_prop)
				);

	// get ready to save PhaseSpace
	if (std::is_same<vfps::interpol_t,float>::value) {
		ps_datatype = H5::PredType::IEEE_F32LE;
	} else {
		ps_datatype = H5::PredType::STD_I32LE;
	}

	static constexpr std::array<hsize_t,ps_rank> ps_maxdims
			= {{ps_xsize,ps_ysize,H5S_UNLIMITED}};

	ps_dataspace = new H5::DataSpace(ps_rank,ps_dims.data(),ps_maxdims.data());


	static constexpr std::array<hsize_t,ps_rank> ps_chunkdims
			= {{ps_xsize/8,ps_ysize/8,1}};
	ps_prop.setChunk(ps_rank,ps_chunkdims.data());
	ps_prop.setShuffle();
	ps_prop.setDeflate(compression);

	ps_dataset = new H5::DataSet(
						file->createDataSet(ps_name,ps_datatype,
											*ps_dataspace,ps_prop)
				);
}

vfps::HDF5File::~HDF5File()
{
	delete file;
	delete bp_dataset;
	delete bp_dataspace;
	delete ps_dataset;
	delete ps_dataspace;
}

void vfps::HDF5File::write(PhaseSpace* ps)
{
	if (std::is_same<vfps::meshdata_t,float>::value) {
		ps_dataset->write(ps->getData(), H5::PredType::NATIVE_FLOAT);
	} else {
		ps_dataset->write(ps->getData(), H5::PredType::NATIVE_INT32);
	}

	if (std::is_same<vfps::integral_t,float>::value) {
		bp_dataset->write(ps->projectionToX(), H5::PredType::NATIVE_FLOAT);
	} else {
		bp_dataset->write(ps->projectionToX(), H5::PredType::NATIVE_INT64);
	}
}

void vfps::HDF5File::append(PhaseSpace* ps)
{
	// append PhaseSpace
	std::array<hsize_t,ps_rank> ps_offset
			= {{0,0,ps_dims[2]}};
	static constexpr std::array<hsize_t,ps_rank> ps_ext
			= {{ps_xsize,ps_ysize,1}};
	ps_dims[2]++;

	ps_dataset->extend(ps_dims.data());
	H5::DataSpace* filespace = new H5::DataSpace(ps_dataset->getSpace());
	filespace->selectHyperslab(H5S_SELECT_SET, ps_ext.data(), ps_offset.data());
	H5::DataSpace* memspace = new H5::DataSpace(ps_rank,ps_ext.data(),nullptr);
	if (std::is_same<vfps::meshdata_t,float>::value) {
		ps_dataset->write(ps->getData(), H5::PredType::NATIVE_FLOAT,
						  *memspace, *filespace);
	} else {
		ps_dataset->write(ps->getData(), H5::PredType::NATIVE_INT32,
						  *memspace, *filespace);
	}

	// append BunchProfile
	std::array<hsize_t,bp_rank> bp_offset
			= {{0,bp_dims[1]}};
	static constexpr std::array<hsize_t,bp_rank> bp_ext
			= {{ps_xsize,1}};
	bp_dims[1]++;

	bp_dataset->extend(bp_dims.data());

	delete filespace;
	filespace = new H5::DataSpace(bp_dataset->getSpace());
	filespace->selectHyperslab(H5S_SELECT_SET, bp_ext.data(), bp_offset.data());

	delete memspace;
	memspace = new H5::DataSpace(bp_rank,bp_ext.data(),nullptr);
	if (std::is_same<vfps::integral_t,float>::value) {
		bp_dataset->write(ps->projectionToX(), H5::PredType::NATIVE_FLOAT,
						  *memspace, *filespace);
	} else {
		bp_dataset->write(ps->projectionToX(), H5::PredType::NATIVE_INT64,
						  *memspace, *filespace);
	}
}
