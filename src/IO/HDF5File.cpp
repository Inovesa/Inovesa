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
#include "IO/HDF5File.hpp"

#include "MessageStrings.hpp"

vfps::HDF5File::HDF5File(const std::string filename,
                         const std::shared_ptr<PhaseSpace> ps,
                         const ElectricField* ef,
                         const std::shared_ptr<Impedance> imp,
                         const WakeFunctionMap* wfm,
                         const size_t nparticles,
                         const double t_sync,
                         const double f_rev)
  : _fname( filename )
  , _file( _prepareFile() )
  , _nBunches( ps->nBunches() )
  , _nParticles( nparticles )
  , _psSizeX( ps->nMeshCells(0) )
  , _psSizeY( ps->nMeshCells(1) )
  , _maxn( (ef != nullptr)? ef->getNMax()/static_cast<size_t>(2) : 0 )
  , _impSize( imp != nullptr ? imp->nFreqs()/2 : 0 )
  , _positionAxis(_makeDatasetInfo<1,meshaxis_t>( "/Info/AxisValues_z"
                                                , {{_psSizeX}}
                                                , {{std::min(2048U,_psSizeX)}}
                                                , {{_psSizeX}}))
  , _energyAxis(_makeDatasetInfo<1,meshaxis_t>( "/Info/AxisValues_E"
                                              , {{_psSizeY}}
                                              , {{std::min(2048U,_psSizeY)}}
                                              , {{_psSizeY}}))
  , _frequencyAxis(_makeDatasetInfo<1,frequency_t>( "/Info/AxisValues_f"
                                                  , {{_maxn}}
                                                  , {{std::min(2048U,_psSizeY)}}
                                                  , {{_maxn}}))
  , _timeAxis(_makeDatasetInfo<1,timeaxis_t>( "/Info/AxisValues_t"
                                             , {{0}},{{256}},{{H5F_UNLIMITED}}))
  , _timeAxisPS(_makeDatasetInfo<1,timeaxis_t>( "/PhaseSpace/axis0"
                                           , {{0}},{{256}},{{H5F_UNLIMITED}}))
  , _bunchPopulation(_makeDatasetInfo<2,integral_t>( "/BunchPopulation/data"
                                                  , {{0,_nBunches}}
                                                  , {{256,_nBunches}}
                                                  , {{H5S_UNLIMITED,_nBunches}}
                                                  ))
  , _bunchProfile(_makeDatasetInfo<3,integral_t>( "/BunchProfile/data"
                                               , {{ 0, _nBunches, _psSizeX }}
                                               , {{ 64, 1
                                                  , std::min(256U,_psSizeX) }}
                                               , {{ H5S_UNLIMITED,_nBunches
                                                  , _psSizeX }} ))
  , _bunchLength( _makeDatasetInfo<2,meshaxis_t>( "/BunchLength/data"
                                               , {{ 0, _nBunches }}
                                               , {{ 256, 1 }}
                                               , {{ H5S_UNLIMITED, _nBunches }}))
  , _bunchPosition( _makeDatasetInfo<2,meshaxis_t>( "/BunchPosition/data"
                                                 , {{ 0, _nBunches }}
                                                 , {{ 256, 1 }}
                                                 , {{ H5S_UNLIMITED, _nBunches }}))
  , _energyProfile(_makeDatasetInfo<3,integral_t>( "/EnergyProfile/data"
                                                 , {{ 0, _nBunches, _psSizeX }}
                                                 , {{ 64, 1
                                                    , std::min(256U,_psSizeX) }}
                                                 , {{ H5S_UNLIMITED,_nBunches
                                                    , _psSizeX }} ))
  , _energySpread( _makeDatasetInfo<2,meshaxis_t>( "/EnergySpread/data"
                                                 , {{ 0, _nBunches }}
                                                 , {{ 256, 1 }}
                                                 , {{ H5S_UNLIMITED
                                                    , _nBunches }}))
  , _energyAverage( _makeDatasetInfo<2,meshaxis_t>( "/EnergyAverage/data"
                                                  , {{ 0, _nBunches }}
                                                  , {{ 256, 1 }}
                                                  , {{ H5S_UNLIMITED, _nBunches }}))
  , _particles( _makeDatasetInfo<3,meshaxis_t>( "/Particles/data"
                                              , {{ 0, _nParticles, 2 }}
                                              , {{ 256
                                                 , std::min(256U,_nParticles)
                                                 , 2 }}
                                              , {{ H5S_UNLIMITED
                                                 , _nParticles, 2 }}))
  , _dynamicRFKick( _makeDatasetInfo<2,meshaxis_t>( "/RFKicks/data"
                                                  , {{ 0, 2 }}
                                                  , {{ 256, 1 }}
                                                  , {{ H5S_UNLIMITED, 2 }}))
  , _wakePotential(_makeDatasetInfo<3,meshaxis_t>( "/WakePotential/data"
                                                 , {{ 0, _nBunches, _psSizeX }}
                                                 , {{ 64, 1
                                                    , std::min(256U,_psSizeX) }}
                                                 , {{ H5S_UNLIMITED,_nBunches
                                                    , _psSizeX }} ))
  , _csrSpectrum(_makeDatasetInfo<3,meshaxis_t>( "/CSR/Spectrum/data"
                                               , {{ 0, _nBunches, _maxn }}
                                               , {{ 64, 1
                                                  , std::min(256U,_maxn) }}
                                               , {{ H5S_UNLIMITED,_nBunches
                                                  , _maxn }} ))
  , _csrIntensity(_makeDatasetInfo<2,meshaxis_t>( "/CSR/Intensity/data"
                                                , {{ 0, _nBunches }}
                                                , {{ 64, 1 }}
                                                , {{ H5S_UNLIMITED
                                                  ,_nBunches }} ))
  , _phaseSpace(_makeDatasetInfo<4,meshdata_t>( "/PhaseSpace/data"
                                              , {{ 0, _nBunches, _psSizeX, _psSizeY }}
                                              , {{ 64, 1
                                                 , std::min(256U,_psSizeX)
                                                 , std::min(256U,_psSizeY) }}
                                              , {{ H5S_UNLIMITED,_nBunches
                                                 , _psSizeX, _psSizeY }} ))
  , _impedanceReal(_makeDatasetInfo<1,csrpower_t>("/Impedance/data/real"
                                                 , {{_impSize}}
                                                 , {{std::min( 4097U,_impSize)}}
                                                 , {{_impSize}}))
  , _impedanceImag(_makeDatasetInfo<1,csrpower_t>("/Impedance/data/imag"
                                                 , {{_impSize}}
                                                 , {{std::min( 4097U,_impSize)}}
                                                 , {{_impSize}}))
  , _wakeFunction(_makeDatasetInfo<1,meshaxis_t>("/WakeFunction/data"
                                                 , {{2*_psSizeX}}
                                                 , {{std::min( 4097U
                                                             , 2U*_psSizeX)}}
                                                 , {{2*_psSizeX}}))
  , _ps(ps)
{
    // Axis
    const double ax_z_meter = ps->getScale(0,"Meter");
    const double ax_z_seconds = ax_z_meter/physcons::c;


    _positionAxis.dataset.createAttribute("Meter",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_z_meter);
    _positionAxis.dataset.createAttribute("Second",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_z_seconds);
    _positionAxis.dataset.write(ps->getAxis(0)->data(),_positionAxis.datatype);

    const double ax_E_eVolt = ps->getScale(1,"ElectronVolt");
    _energyAxis.dataset.createAttribute("ElectronVolt"
                                  , H5::PredType::IEEE_F64LE
                                  , H5::DataSpace()).write(
                                        H5::PredType::IEEE_F64LE,&ax_E_eVolt);
    _energyAxis.dataset.write(ps->getAxis(0)->data(),_energyAxis.datatype);

    if (ef != nullptr || imp != nullptr) {
    // frequency information axis, will be taken from ef or imp
    auto axfreq = (ef != nullptr)? ef->getFreqRuler()
                                 : imp->getRuler();

    const double ax_f_hertz = axfreq->scale("Hertz");
    _frequencyAxis.dataset.createAttribute("Hertz",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_f_hertz);
    _frequencyAxis.dataset.write(axfreq->data(),_frequencyAxis.datatype);
    }

    _timeAxis.dataset.createAttribute("Second",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&t_sync);

    const double axis_t_turns = t_sync*f_rev;
    _timeAxis.dataset.createAttribute("Turn",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&axis_t_turns);

    _bunchPopulation.dataset.createAttribute("Ampere",_bunchPopulation.datatype,
            H5::DataSpace()).write(_bunchPopulation.datatype,&(ps->current));
    _bunchPopulation.dataset.createAttribute("Coulomb",_bunchPopulation.datatype,
            H5::DataSpace()).write(_bunchPopulation.datatype,&(ps->charge));

    // actual data
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/BunchProfile/axis0" );
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/BunchProfile/axis1" );

    _bunchProfile.dataset.createAttribute( "AmperePerNBL"
                                        , H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &ps->current);
    _bunchProfile.dataset.createAttribute( "CoulombPerNBL"
                                        , H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &ps->charge);


    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/BunchLength/axis0" );

    _bunchLength.dataset.createAttribute("Meter",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_z_meter);
    _bunchLength.dataset.createAttribute("Second",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_z_seconds);


    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/BunchPosition/axis0" );

    _bunchPosition.dataset.createAttribute("Meter",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_z_meter);
    _bunchPosition.dataset.createAttribute("Second",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_z_seconds);

    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/EnergyProfile/axis0" );
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_E", "/EnergyProfile/axis1" );

    _energyProfile.dataset.createAttribute("AmperePerNES",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &ps->current);
    _energyProfile.dataset.createAttribute("CoulombPerNES",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &ps->charge);

    _energySpread.dataset.createAttribute("ElectronVolt",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_E_eVolt);

    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/EnergyAverage/axis0" );

    _energyAverage.dataset.createAttribute("ElectronVolt",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_E_eVolt);


    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/Particles/axis0" );


    // get ready to save WakePotential
    if (ef != nullptr) {
        _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/WakePotential/axis0" );
        _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/WakePotential/axis1" );

        _wakePotential.dataset.createAttribute("Volt",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                       &(ef->volts));


        _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/CSR/Spectrum/axis0" );
        _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_f", "/CSR/Spectrum/axis1" );

        _csrIntensity.dataset.createAttribute("WattPerHertz",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                       &ef->factor4WattPerHertz);


        _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/CSR/Intensity/axis0" );

        _csrIntensity.dataset.createAttribute("Watt",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                       &ef->factor4Watts);
    }

    /* The phase space has its own time axis.
     * It resuses information available through
     * variables defined for the original time axis.
     */
    _timeAxisPS.dataset.createAttribute("Second",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&t_sync);
    _timeAxisPS.dataset.createAttribute("Turn",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&axis_t_turns);

    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/PhaseSpace/axis1" );
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_E", "/PhaseSpace/axis2" );

    _phaseSpace.dataset.createAttribute("AmperePerNBLPerNES",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &ps->current);
    _phaseSpace.dataset.createAttribute("CoulombPerNBLPerNES",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &ps->charge);


    if (imp != nullptr) {
        // save Impedance
        _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_f", "/Impedance/axis0" );

        std::vector<csrpower_t> imp_real;
        std::vector<csrpower_t> imp_imag;
        imp_real.reserve(_impSize);
        imp_imag.reserve(_impSize);
        for (const impedance_t& z : imp->impedance()) {
            imp_real.emplace_back(z.real());
            imp_imag.emplace_back(z.imag());
        }
        _impedanceReal.dataset.write(imp_real.data(),_impedanceReal.datatype);
        _impedanceImag.dataset.write(imp_imag.data(),_impedanceImag.datatype);

        _file.openGroup("/Impedance/data").createAttribute("Ohm", H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE, &(imp->factor4Ohms));
    }

    if (wfm != nullptr ) {
        _wakeFunction.dataset.write( wfm->getWakeFunction()
                                   , _wakeFunction.datatype);
    }

    // save Inovesa version
    {
    std::array<hsize_t,1> version_dims {{3}};
    H5::DataSpace version_dspace(1,version_dims.data(),version_dims.data());
    H5::DataSet version_dset = _file.createDataSet
                    ("/Info/Inovesa_v", H5::PredType::STD_I32LE,version_dspace);
    std::array<int32_t,3> version = {{0,0,0}};

    // do not save version for development or feature branches
    auto branchstr = std::string(GIT_BRANCH);
    if (branchstr.empty() || branchstr=="master" || branchstr.front()=='v') {
        version = {{ INOVESA_VERSION_MAJOR,
                     INOVESA_VERSION_MINOR,
                     INOVESA_VERSION_FIX
                  }};
    }
    version_dset.write(version.data(),H5::PredType::NATIVE_INT);
    }

    // save full Inovesa version string
    {
    std::string ver_string = vfps::inovesa_version(true);
    const hsize_t ver_string_dim(ver_string.size());
    H5::DataSpace ver_string_dspace(1,&ver_string_dim);
    H5::DataSet ver_string_dset = _file.createDataSet
                    ("/Info/Inovesa_build", H5::PredType::C_S1,ver_string_dspace);
    ver_string_dset.write(ver_string.c_str(),H5::PredType::C_S1);
    }
}

void vfps::HDF5File::addParameterToGroup(std::string groupname,
                                         std::string paramname,
                                         H5::PredType type,
                                         void* data)
{
    H5::Group group = _file.openGroup(groupname);
    group.createAttribute(paramname,type, H5::DataSpace()).write(type,data);
}

void vfps::HDF5File::append(const ElectricField* ef, const bool fullspectrum)
{
    if (fullspectrum) {
        _appendData(_csrSpectrum,ef->getCSRSpectrum());
    }
    _appendData(_csrIntensity,&ef->getCSRPower());
}

void vfps::HDF5File::appendRFKicks(
        const std::vector<std::array<vfps::meshaxis_t,2>> kicks)
{
    _appendData(_dynamicRFKick,kicks.data(),kicks.size());
}

void vfps::HDF5File::appendTracks(const std::vector<PhaseSpace::Position> &p){
    std::vector<PhaseSpace::Position> physcords;
    physcords.reserve(_nParticles);
    for (auto pos : p) {
        physcords.push_back({_ps->q(pos.x), _ps->p(pos.y)});
    }

    _appendData(_particles,physcords.data());
}

void vfps::HDF5File::append(const PhaseSpace& ps,
                            const timeaxis_t t,
                            const AppendType at)
{
    if ( at == AppendType::All ||
         at == AppendType::PhaseSpace) {
        _appendData(_timeAxisPS,&t);
        _appendData(_phaseSpace,ps.getData());
    }

    if (at != AppendType::PhaseSpace) {
        _appendData(_timeAxis,&t);
        _appendData(_bunchProfile,ps.getProjection(0));
        _appendData(_bunchLength,ps.getBunchLength());
        {
        auto mean_q = ps.getMoment(0,0);
        _appendData(_bunchPosition,mean_q.data());
        }
        _appendData(_energyProfile,ps.getProjection(1));
        _appendData(_energySpread,ps.getEnergySpread());
        {
        auto mean_E = ps.getMoment(1,0);
        _appendData(_energyAverage,mean_E.data());
        }
        {
        _appendData(_bunchPopulation,ps.getBunchPopulation().data());
        }
    }
}

void vfps::HDF5File::append(const vfps::WakeKickMap* wkm)
{
    _appendData(_wakePotential,wkm->getForce());
}



std::unique_ptr<vfps::PhaseSpace>
vfps::HDF5File::readPhaseSpace( std::string fname
                              , meshaxis_t qmin, meshaxis_t qmax
                              , meshaxis_t pmin, meshaxis_t pmax
                              , oclhptr_t oclh
                              , double Qb, double Ib_unscaled
                              , double bl, double dE, int64_t use_step)
{
    H5::DataType datatype;
    if (std::is_same<vfps::meshdata_t,float>::value) {
        datatype = H5::PredType::IEEE_F32LE;
    #if FXP_FRACPART < 31
    } else if (std::is_same<vfps::meshdata_t,fixp32>::value) {
        datatype = H5::PredType::STD_I32LE;
    #endif // FXP_FRACPART < 31
    } else if (std::is_same<vfps::meshdata_t,fixp64>::value) {
        datatype = H5::PredType::STD_I64LE;
    }else if (std::is_same<vfps::meshdata_t,double>::value) {
        datatype = H5::PredType::IEEE_F64LE;
    }

    H5::H5File file(fname,H5F_ACC_RDONLY);

    H5::DataSet ps_dataset = file.openDataSet("/PhaseSpace/data");
    H5::DataSpace ps_space(ps_dataset.getSpace());

    const auto rank = ps_space.getSimpleExtentNdims();
    hsize_t ps_dims[rank];
    ps_space.getSimpleExtentDims( ps_dims, nullptr );

    std::vector<hsize_t> ps_offset;
    std::vector<hsize_t> ps_ext;
    use_step = (ps_dims[0]+use_step)%ps_dims[0];
    meshindex_t ps_size;
    uint32_t nBunches = 1U;
    switch (rank) {
    case 3:
        ps_size = ps_dims[1];
        ps_offset =  {{static_cast<hsize_t>(use_step),0,0}};
        ps_ext = {{1,ps_size,ps_size}};
        break;
    case 4:
        nBunches = ps_dims[1];
        ps_size = ps_dims[2];
        ps_offset =  {{static_cast<hsize_t>(use_step),0,0,0}};
        ps_ext = {{1,nBunches,ps_size,ps_size}};
        break;
    }
    H5::DataSpace memspace(rank,ps_ext.data(),nullptr);
    ps_space.selectHyperslab(H5S_SELECT_SET, ps_ext.data(), ps_offset.data());

    H5::DataType axistype;
    if (std::is_same<vfps::meshaxis_t,float>::value) {
        axistype = H5::PredType::IEEE_F32LE;
    #if FXP_FRACPART < 31
    } else if (std::is_same<vfps::meshaxis_t,fixp32>::value) {
        axistype = H5::PredType::STD_I32LE;
    #endif // FXP_FRACPART < 31
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
        axistype = H5::PredType::STD_I64LE;
    }else if (std::is_same<vfps::meshaxis_t,double>::value) {
        axistype = H5::PredType::IEEE_F64LE;
    }

    auto ps = std::make_unique<PhaseSpace>( ps_size
                                          , qmin,qmax,bl
                                          , pmin,pmax,dE
                                          , oclh
                                          , Qb,Ib_unscaled,1U,1
                                          );
    ps_dataset.read(ps->getData(), datatype, memspace, ps_space);

    return ps;
}

template <int rank, typename datatype>
void vfps::HDF5File::_appendData(DatasetInfo<rank>& ds
                                , const datatype* const data
                                , const size_t size)
{
    std::array<hsize_t,rank> offset{{ds.dims[0]}};
    auto ext = ds.dims;
    ext[0] = size;
    ds.dims[0] += size;
    ds.dataset.extend(ds.dims.data());
    H5::DataSpace filespace(ds.dataset.getSpace());
    filespace.selectHyperslab(H5S_SELECT_SET, ext.data(), offset.data());
    H5::DataSpace memspace(ds.rank,ext.data(),nullptr);
    ds.dataset.write(data, ds.datatype,memspace, filespace);;
}

template <int rank, typename datatype>
vfps::HDF5File::DatasetInfo<rank>
vfps::HDF5File::_makeDatasetInfo( std::string name
                                , std::array<hsize_t,rank> dims
                                , std::array<hsize_t,rank> chunkdims
                                , std::array<hsize_t,rank> maxdims
                                )
{
    for (auto& dim : chunkdims) {
        dim = std::max(dim, static_cast<hsize_t>(1));
    }
    for (auto& dim : maxdims) {
        dim = std::max(dim, static_cast<hsize_t>(1));
    }

    H5::DataType rv_datatype;
    if (std::is_same<datatype,float>::value) {
            rv_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<datatype,double>::value) {
        rv_datatype = H5::PredType::IEEE_F64LE;
    } else if (std::is_same<datatype,fixp64>::value) {
            rv_datatype = H5::PredType::STD_I64LE;
    } else {
        throw std::string("Unknown datatype.");
    }

    H5::DataSpace rv_dataspace(rank,dims.data(),maxdims.data());

    H5::DSetCreatPropList rv_prop;
    rv_prop.setChunk(rank,chunkdims.data());
    rv_prop.setShuffle();
    rv_prop.setDeflate(compression);

    H5::DataSet rv_dataset = _file.createDataSet( name , rv_datatype
                                                , rv_dataspace, rv_prop);


    return DatasetInfo<rank>(rv_dataset,rv_datatype,dims);
}

H5::H5File vfps::HDF5File::_prepareFile()
{
    H5::H5File rv(_fname,H5F_ACC_TRUNC);

    rv.createGroup("/BunchLength");
    rv.createGroup("/BunchPopulation");
    rv.createGroup("/BunchPosition");
    rv.createGroup("/BunchProfile");
    rv.createGroup("/CSR/");
    rv.createGroup("/CSR/Intensity");
    rv.createGroup("/CSR/Spectrum");
    rv.createGroup("/EnergyProfile");
    rv.createGroup("/EnergySpread");
    rv.createGroup("/EnergyAverage");
    rv.createGroup("/Particles");
    rv.createGroup("/PhaseSpace");
    rv.createGroup("/RFKicks");
    rv.createGroup("/Impedance");
    rv.createGroup("/Impedance/data");
    rv.createGroup("/Info");
    rv.createGroup("/Info/Parameters");
    rv.createGroup("/SourceMap");
    rv.createGroup("/SourceMap/data/");
    rv.createGroup("/WakeFunction");
    rv.createGroup("/WakePotential");

    return rv;
}

#endif // INOVESA_USE_HDF5
