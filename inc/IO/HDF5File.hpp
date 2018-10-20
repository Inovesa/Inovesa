// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#pragma once

#if INOVESA_USE_HDF5 == 1

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
     * @brief savePaddedProfile
     * @param ef
     * @param timestep 0 or 1 (for start or end of simulation)
     */
    void appendPadded(const ElectricField* ef);

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

    const uint32_t _nBuckets;

    const uint32_t _nBunches;

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

    DatasetInfo<1> _bucketNumbers;

    DatasetInfo<2> _bunchPopulation;

    DatasetInfo<3> _bunchProfile;

    DatasetInfo<2> _paddedProfile;

    DatasetInfo<2> _bunchLength;

    DatasetInfo<2> _bunchPosition;

    DatasetInfo<3> _energyProfile;

    DatasetInfo<2> _energySpread;

    DatasetInfo<2> _energyAverage;

    DatasetInfo<3> _particles;

    DatasetInfo<2> _dynamicRFKick;

    DatasetInfo<3> _wakePotential;

    DatasetInfo<2> _paddedPotential;

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

#endif // INOVESA_USE_HDF5
