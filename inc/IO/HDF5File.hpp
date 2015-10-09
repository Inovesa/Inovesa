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
#include <gsl/gsl_const_mks.h>
#include <H5Cpp.h>
#include <string>

#include "defines.hpp"
#include "ElectricField.hpp"
#include "PhaseSpace.hpp"

namespace vfps {

class HDF5File
{
public:
	/**
	 * @brief HDF5File
	 * @param fname file name to save HDF5 file to
	 * @param ps_size size of one mesh dimension
	 * @param maxn maximum index (==wavenumber?) of CSR spectrum
	 */
	HDF5File(const std::string fname,
			 const PhaseSpace* ps,
			 const ElectricField* ef);

	~HDF5File();

	void append(const ElectricField* ef);

	void append(const PhaseSpace* ps);

private:
	H5::H5File* file;

	std::string fname;

	static constexpr uint_fast8_t compression = 6;

	H5::DataType datatype_integral;

	H5::DataType datatype_meshdata;

private: // values for phase space axis
	static constexpr uint_fast8_t axps_rank = 2;

	H5::DataSet ax0ps_dataset;

	H5::DataSet ax1ps_dataset;

	H5::DataSpace* ax0ps_dataspace;

	H5::DataSpace* ax1ps_dataspace;

	H5::IntType axps_datatype;

	H5::DSetCreatPropList axps_prop;

private: // values for frequency axis
	static constexpr uint_fast8_t axfreq_rank = 2;

	H5::DataSet axfreq_dataset;

	H5::DataSpace* axfreq_dataspace;

	H5::IntType axfreq_datatype;

	H5::DSetCreatPropList axfreq_prop;


private: // bunch charge
	static constexpr uint_fast8_t bc_rank = 1;

	H5::DataSet bc_dataset;

	H5::DataSpace* bc_dataspace;

	H5::IntType bc_datatype;

	hsize_t bc_dims;

	H5::DSetCreatPropList bc_prop;

private: // bunch profile
	static constexpr uint_fast8_t bp_rank = 2;

	H5::DataSet bp_dataset;

	H5::DataSpace* bp_dataspace;

	H5::IntType bp_datatype;

	std::array<hsize_t,bp_rank> bp_dims;

	H5::DSetCreatPropList bp_prop;

private: // csr spectrum
	static constexpr uint_fast8_t csr_rank = 2;

	H5::DataSet csr_dataset;

	H5::DataSpace* csr_dataspace;

	H5::IntType csr_datatype;

	std::array<hsize_t,csr_rank> csr_dims;

	H5::DSetCreatPropList csr_prop;

	size_t maxn;

private: // phase space
	static constexpr uint_fast8_t ps_rank = 3;

	H5::DataSet ps_dataset;

	H5::DataSpace* ps_dataspace;

	H5::IntType ps_datatype;

	std::array<hsize_t,ps_rank> ps_dims;

	H5::DSetCreatPropList ps_prop;

	meshindex_t ps_size;
};

} // namespace vfps

#endif // HDF5FILE_HPP
