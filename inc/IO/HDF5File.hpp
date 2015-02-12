#ifndef HDF5FILE_HPP
#define HDF5FILE_HPP

#include <array>
#include <H5Cpp.h>
#include <string>

#include "PhaseSpace.hpp"

namespace vfps {

class HDF5File
{
public:
	HDF5File(std::string fname);

	~HDF5File();

	void append(PhaseSpace* ps);

	void write(PhaseSpace* ps);

private:
	H5::H5File* file;

	std::string fname;

	static constexpr unsigned int fs_rank = 3;

	H5::DataSet* fs_dataset;

	H5::DataSpace* fs_dataspace;

	H5::IntType fs_datatype;

	std::array<hsize_t,3> fs_dims;

	const std::string fs_name;

	H5::DSetCreatPropList fs_prop;
};

}

#endif // HDF5FILE_HPP
