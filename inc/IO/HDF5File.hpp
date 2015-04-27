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
	HDF5File(std::string fname, const uint16_t ps_size);

	~HDF5File();

	void append(PhaseSpace* ps);

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

	uint16_t ps_size;
};

}

#endif // HDF5FILE_HPP
