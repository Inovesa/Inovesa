/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
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
#ifndef HDF5FILE_HPP
#define HDF5FILE_HPP

#include <array>
#include <H5Cpp.h>
#include <string>

#include "defines.hpp"
#include "PS/ElectricField.hpp"
#include "Z/Impedance.hpp"
#include "PS/PhaseSpace.hpp"
#include "SM/WakeFunctionMap.hpp"

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
    HDF5File(const std::string filename,
             const PhaseSpace* ps,
             const ElectricField* ef,
             const Impedance* imp,
             const WakeFunctionMap *wfm,
             const double t_sync);

    ~HDF5File();

    void addParameterToGroup(std::string groupname,
                             std::string paramname,
                             H5::PredType type,
                             void* data);

    void append(const ElectricField* ef, bool fullspectrum = true);

    enum class AppendType : uint_fast16_t {
        All, Defaults, PhaseSpace
    };

    void append(const PhaseSpace* ps, const AppendType at=AppendType::Defaults);

    void append(const WakeKickMap* wkm);

    void appendTime(const double t);

public:
    static PhaseSpace readPhaseSpace(std::string fname,
                                     meshaxis_t qmin, meshaxis_t qmax,
                                     meshaxis_t pmin, meshaxis_t pmax,
                                     double Qb, double Ib_unscaled,
                                     double bl, double dE,
                                     int64_t use_step=-1l);

private:
    H5::H5File* _file;

    std::string fname;

    static constexpr uint_fast8_t compression = 6;

    H5::DataType datatype_integral;

    H5::DataType datatype_meshdata;

private: // values for phase space axis
    static constexpr uint_fast8_t axps_rank = 1;

    H5::DataSet ax0ps_dataset;

    H5::DataSet ax1ps_dataset;

    H5::DataType axps_datatype;

    H5::DSetCreatPropList axps_prop;

private: // values for frequency axis
    static constexpr uint_fast8_t axfreq_rank = 1;

    H5::DataSet axfreq_dataset;

    H5::DataType axfreq_datatype;

    H5::DSetCreatPropList axfreq_prop;

private: // time axis
    static constexpr uint_fast8_t ta_rank = 1;

    H5::DataSet ta_dataset;

    H5::DataType ta_datatype;

    hsize_t ta_dims;

    H5::DSetCreatPropList ta_prop;

private: // bunch current
    static constexpr uint_fast8_t bc_rank = 1;

    H5::DataSet bc_dataset;

    H5::DataType bc_datatype;

    hsize_t bc_dims;

    H5::DSetCreatPropList bc_prop;

private: // bunch profile
    static constexpr uint_fast8_t bp_rank = 2;

    H5::DataSet bp_dataset;

    H5::DataType bp_datatype;

    std::array<hsize_t,bp_rank> bp_dims;

    H5::DSetCreatPropList bp_prop;

private: // bunch length
    static constexpr uint_fast8_t bl_rank = 1;

    H5::DataSet bl_dataset;

    H5::DataType bl_datatype;

    hsize_t bl_dims;

    H5::DSetCreatPropList bl_prop;

private: // bunch position
    static constexpr uint_fast8_t qb_rank = 1;

    H5::DataSet qb_dataset;

    H5::DataType qb_datatype;

    hsize_t qb_dims;

    H5::DSetCreatPropList qb_prop;

private: // energy profile
    static constexpr uint_fast8_t ep_rank = 2;

    H5::DataSet ep_dataset;

    H5::DataType ep_datatype;

    std::array<hsize_t,ep_rank> ep_dims;

    H5::DSetCreatPropList ep_prop;

private: // energy spread
    static constexpr uint_fast8_t es_rank = 1;

    H5::DataSet es_dataset;

    H5::DataType es_datatype;

    hsize_t es_dims;

    H5::DSetCreatPropList es_prop;

private: // wake potential
    static constexpr uint_fast8_t wp_rank = 2;

    H5::DataSet wp_dataset;

    H5::DataType wp_datatype;

    std::array<hsize_t,wp_rank> wp_dims;

    H5::DSetCreatPropList wp_prop;

private: // csr spectrum
    static constexpr uint_fast8_t csr_rank = 2;

    H5::DataSet csr_dataset;

    H5::DataType csr_datatype;

    std::array<hsize_t,csr_rank> csr_dims;

    H5::DSetCreatPropList csr_prop;

    size_t maxn;

private: // csr intensity
    static constexpr uint_fast8_t csri_rank = 1;

    H5::DataSet csri_dataset;

    H5::DataType csri_datatype;

    hsize_t csri_dims;

    H5::DSetCreatPropList csri_prop;

private: // phase space
    static constexpr uint_fast8_t ps_rank = 3;

    H5::DataSet _ps_dataset;

    H5::DataType ps_datatype;

    std::array<hsize_t,ps_rank> ps_dims;

    H5::DSetCreatPropList ps_prop;

    meshindex_t ps_size;

private: // impedance
    static constexpr uint_fast8_t imp_rank = 1;

    H5::DataSet imp_dataset_real;
    H5::DataSet imp_dataset_imag;

    H5::DataType imp_datatype;

    hsize_t imp_size;

    H5::DSetCreatPropList imp_prop;

private: // wake function
    static constexpr uint_fast8_t wf_rank = 1;

    H5::DataSet wf_dataset;

    H5::DataType wf_datatype;

    hsize_t wf_size;

    H5::DSetCreatPropList wf_prop;
};

} // namespace vfps

#endif // HDF5FILE_HPP
#endif // INOVESA_USE_HDF5
