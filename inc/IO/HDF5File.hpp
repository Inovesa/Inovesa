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

    void appendRFKicks(const std::vector<std::array<vfps::meshaxis_t,2>> kicks);

    void appendTracks(const std::vector<PhaseSpace::Position> &p);

    void append(const PhaseSpace& ps,
                const timeaxis_t t,
                const AppendType at=AppendType::Defaults);

    void append(const WakeKickMap* wkm);

public:
    static std::unique_ptr<PhaseSpace>
        readPhaseSpace( std::string fname
                      , meshaxis_t qmin, meshaxis_t qmax
                      , meshaxis_t pmin, meshaxis_t pmax
                      , oclhptr_t oclh
                      , double Qb, double Ib_unscaled
                      , double bl, double dE
                      , int64_t use_step=-1l
                      );

private:
    template <int N>
    struct DatasetInfo {
        DatasetInfo() = delete;
        DatasetInfo( H5::DataSet dataset
                   , H5::DataType datatype
                   , std::array<hsize_t,N> dims)
        : dataset(dataset), datatype(datatype), dims(dims)
        {}

        static constexpr int rank = N;
        const H5::DataSet dataset;
        const H5::DataType datatype;
        std::array<hsize_t,N> dims;
    };

    const std::string _fname;

    H5::H5File _file;

    const meshindex_t _nBunches;

    const uint32_t _nParticles;

    const uint32_t _psSizeX;

    const uint32_t _psSizeY;

    const uint32_t _maxn;

    const uint32_t _impSize;

    static constexpr uint_fast8_t compression = 6;

private: // values for phase space axis

    DatasetInfo<1> _positionAxis;

    DatasetInfo<1> _energyAxis;

    DatasetInfo<1> _frequencyAxis;

    DatasetInfo<1> _timeAxis;

    DatasetInfo<1> _timeAxisPS;

    DatasetInfo<2> _bunchPopulation;

    DatasetInfo<3> _bunchProfile;

    DatasetInfo<2> _bunchLength;

    DatasetInfo<2> _bunchPosition;

    DatasetInfo<3> _energyProfile;

    DatasetInfo<2> _energySpread;

    DatasetInfo<2> _energyAverage;

    DatasetInfo<3> _particles;

    DatasetInfo<2> _dynamicRFKick;

    DatasetInfo<3> _wakePotential;

    DatasetInfo<3> _csrSpectrum;

    DatasetInfo<2> _csrIntensity;

    DatasetInfo<4> _phaseSpace;

    DatasetInfo<1> _impedanceReal;
    DatasetInfo<1> _impedanceImag;

    DatasetInfo<1> _wakeFunction;

private:
    template <int rank, typename datatype>
    void _appendData( DatasetInfo<rank>& ds
                    , const datatype* const data
                    , const size_t size=1);

    template<int rank, typename datatype>
    DatasetInfo<rank> _makeDatasetInfo( std::string name
                                      , std::array<hsize_t,rank> dims
                                      , std::array<hsize_t,rank> chunkdims
                                      , std::array<hsize_t,rank> maxdims);


    H5::H5File _prepareFile();

private:
    std::shared_ptr<PhaseSpace> _ps;
};

} // namespace vfps

#endif // HDF5FILE_HPP
#endif // INOVESA_USE_HDF5
