#include "IO/HDF5File.hpp"

vfps::HDF5File::HDF5File(std::string fname) :
	file( nullptr ),
	fname( fname ),
	fs_datatype( H5::PredType::IEEE_F32LE ),
	fs_dims( {{ps_xsize,ps_ysize,1}} ),
	fs_name( "PhaseSpace" )
{
	file = new H5::H5File(fname,H5F_ACC_TRUNC);

	static constexpr std::array<hsize_t,fs_rank> fs_maxdims
			= {{ps_xsize,ps_ysize,H5S_UNLIMITED}};

	fs_dataspace = new H5::DataSpace(fs_rank,fs_dims.data(),fs_maxdims.data());


	static constexpr std::array<hsize_t,fs_rank> fs_chunkdims
			= {{ps_xsize/8,ps_ysize/8,1}};
	fs_prop.setChunk(3,fs_chunkdims.data());
	fs_prop.setShuffle();
	fs_prop.setDeflate(6);

	fs_dataset = new H5::DataSet(
						file->createDataSet(fs_name,fs_datatype,
											*fs_dataspace,fs_prop)
				);
}

vfps::HDF5File::~HDF5File()
{
	delete file;
	delete fs_dataset;
	delete fs_dataspace;
}

void vfps::HDF5File::write(PhaseSpace* ps)
{
	fs_dataset->write(ps->getData(), H5::PredType::NATIVE_FLOAT);
}

void vfps::HDF5File::append(PhaseSpace* ps)
{
	std::array<hsize_t,fs_rank> offset
			= {{0,0,fs_dims[2]}};
	static constexpr std::array<hsize_t,fs_rank> fs_ext
			= {{ps_xsize,ps_ysize,1}};
	fs_dims[2]++;

	fs_dataset->extend(fs_dims.data());
	H5::DataSpace* filespace = new H5::DataSpace(fs_dataset->getSpace());
	filespace->selectHyperslab(H5S_SELECT_SET, fs_ext.data(), offset.data());
	H5::DataSpace* memspace = new H5::DataSpace(fs_rank,fs_ext.data(),nullptr);
	fs_dataset->write(ps->getData(), H5::PredType::NATIVE_FLOAT,
					  *memspace, *filespace);

	delete memspace;
	delete filespace;
}
