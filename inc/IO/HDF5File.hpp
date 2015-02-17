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

	static constexpr unsigned int compression = 6;

private:
	static constexpr unsigned int bp_rank = 2;

	H5::DataSet* bp_dataset;

	H5::DataSpace* bp_dataspace;

	H5::IntType bp_datatype;

	std::array<hsize_t,bp_rank> bp_dims;

	const std::string bp_name;

	H5::DSetCreatPropList bp_prop;

private:
	static constexpr unsigned int ps_rank = 3;

	H5::DataSet* ps_dataset;

	H5::DataSpace* ps_dataspace;

	H5::IntType ps_datatype;

	std::array<hsize_t,ps_rank> ps_dims;

	const std::string ps_name;

	H5::DSetCreatPropList ps_prop;
};

}

#endif // HDF5FILE_HPP
