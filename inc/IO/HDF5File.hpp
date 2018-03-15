/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
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
#include <memory>
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
     * @param filename file name to save HDF5 file to
     * @param ps phase space (to know dimensions, etc.)
     * @param ef electric field used for CSR computation
     * @param imp Impedance
     * @param wfm wake
     * @param nparticles
     * @param t_sync
     */
    HDF5File(const std::string filename,
             const std::shared_ptr<PhaseSpace> ps,
             const ElectricField* ef,
             const std::shared_ptr<Impedance> imp,
             const WakeFunctionMap *wfm,
             const size_t nparticles,
             const double t_sync,
             const double f_rev);

    ~HDF5File() = default;

    void addParameterToGroup(std::string groupname,
                             std::string paramname,
                             H5::PredType type,
                             void* data);

    /**
     * @brief append synchrotron radiation data
     * @param ef electric fild used for CSR computation
     * @param fullspectrum safe full CSR spectrum or only integrated power
     */
    void append(const ElectricField* ef, const bool fullspectrum = true);

    /**
     * @brief The AppendType enum
     *
     * All: save everything
     * Defaults: (no phase space)
     * PhaseSpace: phase space only
     */
    enum class AppendType : uint_fast16_t {
        All, Defaults, PhaseSpace
    };

    void appendTracks(const PhaseSpace::Position *particles);

    void appendSourceMap(const PhaseSpace::Position *allpos);

    void append(const PhaseSpace& ps,
                const timeaxis_t t,
                const AppendType at=AppendType::Defaults);

    void append(const WakeKickMap* wkm);

public:
    static std::unique_ptr<PhaseSpace>
        readPhaseSpace(std::string fname
                      , meshaxis_t qmin, meshaxis_t qmax
                      , meshaxis_t pmin, meshaxis_t pmax
                      #ifdef INOVESA_USE_OPENCL
                      , std::shared_ptr<OCLH> oclh
                      #endif // INOVESA_USE_OPENCL
                      , double Qb, double Ib_unscaled
                      , double bl, double dE
                      , int64_t use_step=-1l
                      );

private:
    std::unique_ptr<H5::H5File> _file;

    std::string fname;

    static constexpr uint_fast8_t compression = 6;

    H5::DataType datatype_integral;

    H5::DataType datatype_meshdata;

private: // values for phase space axis
    static constexpr uint_fast8_t axps_rank = 1;

    H5::DataSet ax0ps_dataset;

    H5::DataSet ax1ps_dataset;

    H5::DataType axps_datatype;

private: // values for frequency axis
    static constexpr uint_fast8_t axfreq_rank = 1;

    H5::DataSet axfreq_dataset;

    H5::DataType axfreq_datatype;

private: // time axis
    static constexpr uint_fast8_t ta_rank = 1;

    H5::DataSet ta_dataset;

    H5::DataType ta_datatype;

    hsize_t ta_dims;

private: // bunch current
    static constexpr uint_fast8_t bc_rank = 1;

    H5::DataSet bc_dataset;

    H5::DataType bc_datatype;

    hsize_t bc_dims;

private: // bunch profile
    static constexpr uint_fast8_t bp_rank = 2;

    H5::DataSet bp_dataset;

    H5::DataType bp_datatype;

    std::array<hsize_t,bp_rank> bp_dims;

private: // bunch length
    static constexpr uint_fast8_t bl_rank = 1;

    H5::DataSet bl_dataset;

    H5::DataType bl_datatype;

    hsize_t bl_dims;

private: // bunch position
    static constexpr uint_fast8_t qb_rank = 1;

    H5::DataSet qb_dataset;

    H5::DataType qb_datatype;

    hsize_t qb_dims;

private: // energy profile
    static constexpr uint_fast8_t ep_rank = 2;

    H5::DataSet ep_dataset;

    H5::DataType ep_datatype;

    std::array<hsize_t,ep_rank> ep_dims;

private: // energy spread
    static constexpr uint_fast8_t es_rank = 1;

    H5::DataSet es_dataset;

    H5::DataType es_datatype;

    hsize_t es_dims;

private: // particles (tracking)
    static constexpr uint_fast8_t pt_rank = 3;

    H5::DataSet pt_dataset;

    H5::DataType pt_datatype;

    std::array<hsize_t,pt_rank> pt_dims;

    const size_t pt_particles;

private: // wake potential
    static constexpr uint_fast8_t wp_rank = 2;

    H5::DataSet wp_dataset;

    H5::DataType wp_datatype;

    std::array<hsize_t,wp_rank> wp_dims;

private: // csr spectrum
    static constexpr uint_fast8_t csr_rank = 2;

    H5::DataSet csr_dataset;

    H5::DataType csr_datatype;

    size_t maxn;

    std::array<hsize_t,csr_rank> csr_dims;

private: // csr intensity
    static constexpr uint_fast8_t csri_rank = 1;

    H5::DataSet csri_dataset;

    H5::DataType csri_datatype;

    hsize_t csri_dims;

private: // phase space
    static constexpr uint_fast8_t _ps_rank = 3;

    H5::DataSet _ps_dataset;

    H5::DataSet _ps_ta_dataset;

    H5::DataType _ps_datatype;

    std::array<hsize_t,_ps_rank> _ps_dims;

    meshindex_t _ps_size;

    hsize_t _ps_ta_dims;

private: // impedance
    static constexpr uint_fast8_t imp_rank = 1;

    H5::DataSet imp_dataset_real;
    H5::DataSet imp_dataset_imag;

    H5::DataType imp_datatype;

    hsize_t imp_size;

private: // source map
    static constexpr uint_fast8_t _sm_rank = _ps_rank;

    H5::DataSet _sm_dataset_x;
    H5::DataSet _sm_dataset_y;

    H5::DataType _sm_datatype;

    std::array<hsize_t,_sm_rank> _sm_dims;

    meshindex_t _sm_size;

private: // wake function
    static constexpr uint_fast8_t wf_rank = 1;

    H5::DataSet wf_dataset;

    H5::DataType wf_datatype;

    hsize_t wf_size;
};

} // namespace vfps

#endif // HDF5FILE_HPP
#endif // INOVESA_USE_HDF5
