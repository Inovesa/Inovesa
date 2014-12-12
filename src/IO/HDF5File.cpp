#include "IO/HDF5File.hpp"

vfps::HDF5File::HDF5File(std::string fname) :
	file( nullptr ),
	fname( fname ),
	fs_datatype( H5::PredType::IEEE_F32LE ),
	fs_dims( {fs_xsize,fs_ysize,1} ),
	fs_name( "PhaseSpace" )
{
	file = new H5::H5File(fname,H5F_ACC_TRUNC);

	static constexpr std::array<hsize_t,fs_rank> fs_maxdims
			= {fs_xsize,fs_ysize,H5S_UNLIMITED};

	fs_dataspace = new H5::DataSpace(fs_rank,fs_dims.data(),fs_maxdims.data());


	static constexpr std::array<hsize_t,fs_rank> fs_chunkdims
			= {fs_xsize/8,fs_ysize/8,1};
	fs_prop.setChunk(3,fs_chunkdims.data());

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

void vfps::HDF5File::write(Mesh2D<meshdata_t>* ps)
{
	fs_dataset->write(ps->getData(), H5::PredType::NATIVE_FLOAT);
}

void vfps::HDF5File::append(Mesh2D<meshdata_t>* ps)
{
	static constexpr std::array<hsize_t,fs_rank> offset
			= {0,0,1};
	static constexpr std::array<hsize_t,fs_rank> fs_ext
			= {fs_xsize,fs_ysize,1};
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
