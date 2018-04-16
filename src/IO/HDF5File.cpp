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
  , _psSizeX( ps->nMeshCells(0) )
  , _psSizeY( ps->nMeshCells(1) )
  , timeaxis(_makeDatasetInfo<1,timeaxis_t>( "/Info/AxisValues_t"
                                           , {{0}},{{256}},{{H5F_UNLIMITED}}))
  , bunchpopulation(_makeDatasetInfo<2,integral_t>( "/BunchPopulation/data"
                                                  , {{0,_nBunches}}
                                                  , {{256,1}}
                                                  , {{H5S_UNLIMITED,_nBunches}}
                                                  ))
  , bunchprofile(_makeDatasetInfo<3,integral_t>( "/BunchProfile/data"
                                               , {{ 0, _nBunches, _psSizeX }}
                                               , {{ 64, 1
                                                  , std::min(256U,_psSizeX) }}
                                               , {{ H5S_UNLIMITED,_nBunches
                                                  , _psSizeX }} ))
  , bunchlength( _makeDatasetInfo<2,meshaxis_t>( "/BunchLength/data"
                                               , {{ 0, _nBunches }}
                                               , {{ 256, 1 }}
                                               , {{ H5S_UNLIMITED, _nBunches }}))
  , bunchposition( _makeDatasetInfo<2,meshaxis_t>( "/BunchPosition/data"
                                                 , {{ 0, _nBunches }}
                                                 , {{ 256, 1 }}
                                                 , {{ H5S_UNLIMITED, _nBunches }}))
  , ep_dims( {{ 0, ps->nBunches(), ps->nMeshCells(0) }} )
  , es_dims( {{ 0, ps->nBunches() }} )
  , pt_dims( {{ 0, nparticles, 2 }} )
  , pt_particles( nparticles )
  , drfk_dims({{ 0,2 }} )
  , wp_dims( {{ 0, ps->nBunches(), ps->nMeshCells(0) }} )
  , maxn( (ef != nullptr)? ef->getNMax()/static_cast<size_t>(2) : 0 )
  , csr_dims( {{ 0, ps->nBunches(), maxn }} )
  , csri_dims( {{ 0, ps->nBunches() }} )
  , _ps_dims( {{ 0, ps->nBunches(), ps->nMeshCells(0), ps->nMeshCells(1) }} )
  , _ps_ta_dims( 0 )
  , imp_size( imp != nullptr ? imp->nFreqs()/2 : 0 )
  , _sm_dims (_ps_dims )
  , wf_size( 2*_psSizeX )
{
    // save Values of Phase Space Axis
    if (std::is_same<vfps::meshaxis_t,float>::value) {
            axps_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            axps_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            axps_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,axps_rank> psa_dims
                    = {{ _psSizeX }};
    const std::array<hsize_t,axps_rank> psa_maxdims
                    = {{ _psSizeX }};

    H5::DataSpace ax0ps_dataspace(axps_rank,psa_dims.data(),
                                        psa_maxdims.data());
    H5::DataSpace ax1ps_dataspace(axps_rank,psa_dims.data(),
                                        psa_maxdims.data());


    const hsize_t psa_chunkdims = std::min(2048U,_psSizeX);

    H5::DSetCreatPropList axps_prop;
    axps_prop.setChunk(axps_rank,&psa_chunkdims);
    axps_prop.setShuffle();
    axps_prop.setDeflate(compression);

    ax0ps_dataset = _file.createDataSet("/Info/AxisValues_z",axps_datatype,
                                         ax0ps_dataspace,axps_prop);
    const double ax_z_meter = ps->getScale(0);
    ax0ps_dataset.createAttribute("Meter",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_z_meter);
    const double ax_z_seconds = ax_z_meter/physcons::c;
    ax0ps_dataset.createAttribute("Second",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_z_seconds);

    ax1ps_dataset = _file.createDataSet("/Info/AxisValues_E",axps_datatype,
                                         ax1ps_dataspace,axps_prop);
    const double ax1scale = ps->getScale(1);
    ax1ps_dataset.createAttribute("ElectronVolt"
                                  ,H5::PredType::IEEE_F64LE,
                                  H5::DataSpace()).write(
          H5::PredType::IEEE_F64LE,&ax1scale);
    ax0ps_dataset.write(ps->getAxis(0)->data(),axps_datatype);
    ax1ps_dataset.write(ps->getAxis(0)->data(),axps_datatype);

    // frequency information axis, will be taken from ef or imp
    {
    const Ruler<frequency_t>* axfreq{nullptr};

    if (ef != nullptr) {
        axfreq = ef->getFreqRuler();
    } else if (imp != nullptr) {
        axfreq = imp->getRuler();
    }

    // save Values of Frequency Axis
    if (std::is_same<vfps::meshaxis_t,float>::value) {
            axfreq_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            axfreq_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            axfreq_datatype = H5::PredType::IEEE_F64LE;
    }
    const hsize_t axfreq_dims = maxn;
    const hsize_t axfreq_maxdims =  maxn;

    H5::DataSpace axfreq_dataspace(axfreq_rank,&axfreq_dims,&axfreq_maxdims);

    const hsize_t axfreq_chunkdims = std::min(static_cast<size_t>(2048),
                                              maxn);


    H5::DSetCreatPropList axfreq_prop;
    axfreq_prop.setChunk(axps_rank,&axfreq_chunkdims);
    axfreq_prop.setShuffle();
    axfreq_prop.setDeflate(compression);

    axfreq_dataset = _file.createDataSet("/Info/AxisValues_f",axfreq_datatype,
                                          axfreq_dataspace,axfreq_prop);
    const double axfreqscale = axfreq->scale();
    axfreq_dataset.createAttribute("Hertz",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&axfreqscale);
    axfreq_dataset.write(axfreq->data(),axfreq_datatype);
    }

    timeaxis.dataset.createAttribute("Second",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&t_sync);

    const double axturnscale = t_sync*f_rev;
    timeaxis.dataset.createAttribute("Turn",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&axturnscale);

    bunchpopulation.dataset.createAttribute("Ampere",bunchpopulation.datatype,
            H5::DataSpace()).write(bunchpopulation.datatype,&(ps->current));
    bunchpopulation.dataset.createAttribute("Coulomb",bunchpopulation.datatype,
            H5::DataSpace()).write(bunchpopulation.datatype,&(ps->charge));


    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/BunchProfile/axis0" );
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/BunchProfile/axis1" );

    bunchprofile.dataset.createAttribute( "AmperePerNBL"
                                        , H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &ps->current);
    bunchprofile.dataset.createAttribute( "CoulombPerNBL"
                                        , H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &ps->charge);


    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/BunchLength/axis0" );

    bunchlength.dataset.createAttribute("Meter",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_z_meter);
    bunchlength.dataset.createAttribute("Second",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_z_seconds);


    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/BunchPosition/axis0" );

    bunchposition.dataset.createAttribute("Meter",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_z_meter);
    bunchposition.dataset.createAttribute("Second",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax_z_seconds);


    // get ready to save EnergyProfiles
    {
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/EnergyProfile/axis0" );
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_E", "/EnergyProfile/axis1" );

    if (std::is_same<vfps::integral_t,float>::value) {
            ep_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::integral_t,fixp64>::value) {
            ep_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::integral_t,double>::value) {
            ep_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,ep_rank> ep_maxdims
                    = {{H5S_UNLIMITED,_nBunches, _psSizeX}};

    H5::DataSpace ep_dataspace(ep_rank,ep_dims.data(),ep_maxdims.data());

    const std::array<hsize_t,ep_rank> ep_chunkdims
        = {{64U,1U,std::min(2048U,_psSizeX)}};

    H5::DSetCreatPropList ep_prop;
    ep_prop.setChunk(ep_rank,ep_chunkdims.data());
    ep_prop.setShuffle();
    ep_prop.setDeflate(compression);


    const double ep_AmperePerSigma = ps->current;
    const double ep_CoulombPerSigma = ps->charge;

    ep_dataset = _file.createDataSet("/EnergyProfile/data",ep_datatype,
                                      ep_dataspace,ep_prop);

    ep_dataset.createAttribute("AmperePerNES",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &ep_AmperePerSigma);
    ep_dataset.createAttribute("CoulombPerNES",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &ep_CoulombPerSigma);
    }

    // get ready to save Energy Spread
    {
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/EnergySpread/axis0" );
    if (std::is_same<vfps::meshaxis_t,float>::value) {
            es_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            es_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            es_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,es_rank> es_maxdims {{H5S_UNLIMITED, _nBunches }};

    H5::DataSpace es_dataspace(es_rank,es_dims.data(),es_maxdims.data());

    const std::array<hsize_t,es_rank> es_chunkdims {{ 256U, 1U }};

    H5::DSetCreatPropList es_prop;
    es_prop.setChunk(es_rank,es_chunkdims.data());
    es_prop.setShuffle();
    es_prop.setDeflate(compression);

    es_dataset = _file.createDataSet("/EnergySpread/data",es_datatype,
                                      es_dataspace,es_prop);
    es_dataset.createAttribute("ElectronVolt",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax1scale);
    }

    // get ready to save particles from (pseudo-) tracking
    {
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/Particles/axis0" );
    if (std::is_same<vfps::meshaxis_t,float>::value) {
            pt_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            pt_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            pt_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,pt_rank> pt_maxdims
            = {{H5S_UNLIMITED,
                std::max(nparticles,static_cast<decltype(nparticles)>(1U)),
                2U}};

    H5::DataSpace pt_dataspace(pt_rank,pt_dims.data(),pt_maxdims.data());

    const std::array<hsize_t,pt_rank> pt_chunkdims
            = {{64U,
                std::max(static_cast<decltype(nparticles)>(1U),
                std::min(static_cast<decltype(nparticles)>(1024U),nparticles)),
                2U}};

    H5::DSetCreatPropList pt_prop;
    pt_prop.setChunk(pt_rank,pt_chunkdims.data());
    pt_prop.setShuffle();
    pt_prop.setDeflate(compression);

    pt_dataset = _file.createDataSet("/Particles/data",pt_datatype,
                                      pt_dataspace,pt_prop);
    }

    // get ready to save status of dynamic RF kick
    {
    if (std::is_same<vfps::meshaxis_t,float>::value) {
            drfk_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            drfk_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            drfk_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,drfk_rank> drfk_maxdims = {{H5S_UNLIMITED,2U}};

    H5::DataSpace drfk_dataspace(drfk_rank,drfk_dims.data(),drfk_maxdims.data());

    const std::array<hsize_t,drfk_rank> drfk_chunkdims = {{64U,2U}};

    H5::DSetCreatPropList drfk_prop;
    drfk_prop.setChunk(drfk_rank,drfk_chunkdims.data());
    drfk_prop.setShuffle();
    drfk_prop.setDeflate(compression);

    drfk_dataset = _file.createDataSet( "/RFKicks/data",drfk_datatype
                                       , drfk_dataspace,drfk_prop);
    }

    // get ready to save WakePotential
    if (ef != nullptr) {
        _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/WakePotential/axis0" );
        _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/WakePotential/axis1" );

        if (std::is_same<vfps::meshaxis_t,float>::value) {
            wp_datatype = H5::PredType::IEEE_F32LE;
        } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            wp_datatype = H5::PredType::STD_I64LE;
        } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            wp_datatype = H5::PredType::IEEE_F64LE;
        }

        const std::array<hsize_t,wp_rank> wp_maxdims
                = {{H5S_UNLIMITED,_nBunches,_psSizeX}};

        H5::DataSpace wp_dataspace(wp_rank,wp_dims.data(),wp_maxdims.data());

        const std::array<hsize_t,wp_rank> wp_chunkdims
                = {{64U,1U,std::min(2048U,_psSizeX)}};

        H5::DSetCreatPropList wp_prop;
        wp_prop.setChunk(wp_rank,wp_chunkdims.data());
        wp_prop.setShuffle();
        wp_prop.setDeflate(compression);

        wp_dataset = _file.createDataSet("/WakePotential/data",wp_datatype,
                                          wp_dataspace,wp_prop);
        wp_dataset.createAttribute("Volt",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                       &(ef->volts));

        // get ready to save CSR Data
        // get ready to save CSR Spectrum
        _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/CSR/Spectrum/axis0" );
        _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_f", "/CSR/Spectrum/axis1" );

        if (std::is_same<vfps::csrpower_t,float>::value) {
            csr_datatype = H5::PredType::IEEE_F32LE;
        } else if (std::is_same<vfps::csrpower_t,double>::value) {
            csr_datatype = H5::PredType::IEEE_F64LE;
        }

        const std::array<hsize_t,csr_rank> csr_maxdims
            = {{H5S_UNLIMITED,_nBunches,maxn}};

        H5::DataSpace csr_dataspace(csr_rank,csr_dims.data(),csr_maxdims.data());

        const std::array<hsize_t,csr_rank> csr_chunkdims
            = {{64U,1U,std::min(static_cast<size_t>(2048),maxn)}};

        H5::DSetCreatPropList csr_prop;
        csr_prop.setChunk(csr_rank,csr_chunkdims.data());
        csr_prop.setShuffle();
        csr_prop.setDeflate(compression);

        const double csrs_WattPerHertz = ef->factor4WattPerHertz;

        csr_dataset = _file.createDataSet("/CSR/Spectrum/data",csr_datatype,
                                           csr_dataspace,csr_prop);
        csr_dataset.createAttribute("WattPerHertz",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                       &csrs_WattPerHertz);

        // get ready to save CSR Intensity
        _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/CSR/Intensity/axis0" );

        if (std::is_same<vfps::csrpower_t,float>::value) {
            csri_datatype = H5::PredType::IEEE_F32LE;
        } else if (std::is_same<vfps::csrpower_t,double>::value) {
            csri_datatype = H5::PredType::IEEE_F64LE;
        }

        const std::array<hsize_t,csri_rank> csri_maxdims {{H5S_UNLIMITED
                                                          , _nBunches}};

        H5::DataSpace csrp_dataspace( csri_rank,csri_dims.data()
                                    , csri_maxdims.data());


        const std::array<hsize_t,csri_rank>
                csri_chunkdims{{std::min(2048U,_psSizeX),1U}};

        H5::DSetCreatPropList csri_prop;
        csri_prop.setChunk(csri_rank,csri_chunkdims.data());
        csri_prop.setShuffle();
        csri_prop.setDeflate(compression);

        csri_dataset = _file.createDataSet("/CSR/Intensity/data",csri_datatype,
                                            csrp_dataspace,csri_prop);

        const double csri_watt = ef->factor4Watts;

        csri_dataset.createAttribute("Watt",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                       &csri_watt);
    }

    // get ready to save PhaseSpace
    {

    /* The phase space has its own time axis.
     * It resuses information available through
     * variables defined in the headerfrom the original time axis.
     */

    hsize_t ta_maxdims = H5S_UNLIMITED;

    H5::DataSpace ta_dataspace(timeaxis.rank,&_ps_ta_dims,&ta_maxdims);


    const hsize_t ta_chunkdims = 256;

    H5::DSetCreatPropList ta_prop;
    ta_prop.setChunk(timeaxis.rank,&ta_chunkdims);
    ta_prop.setShuffle();
    ta_prop.setDeflate(compression);

    _ps_ta_dataset = _file.createDataSet("/PhaseSpace/axis0",timeaxis.datatype,
                                          ta_dataspace,ta_prop);
    const double axtimescale = t_sync;
    _ps_ta_dataset.createAttribute("Second",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&axtimescale);

    const double axturnscale = axtimescale*f_rev;
    _ps_ta_dataset.createAttribute("Turn",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&axturnscale);

    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/PhaseSpace/axis1" );
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_E", "/PhaseSpace/axis2" );

    if (std::is_same<vfps::meshdata_t,float>::value) {
        _ps_datatype = H5::PredType::IEEE_F32LE;
    #if FXP_FRACPART < 31
    } else if (std::is_same<vfps::meshdata_t,fixp32>::value) {
        _ps_datatype = H5::PredType::STD_I32LE;
    #endif // FXP_FRACPART < 31
    } else if (std::is_same<vfps::meshdata_t,fixp64>::value) {
        _ps_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshdata_t,double>::value) {
        _ps_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,_ps_rank> ps_maxdims
            = {{H5S_UNLIMITED,_nBunches,_psSizeX,_psSizeY}};

    H5::DataSpace ps_dataspace(_ps_rank,_ps_dims.data(),ps_maxdims.data());


    const std::array<hsize_t,_ps_rank> ps_chunkdims
        = {{64U,1U,std::min(64U,_psSizeX),std::min(64U,_psSizeX)}};

    H5::DSetCreatPropList ps_prop;
    ps_prop.setChunk(_ps_rank,ps_chunkdims.data());
    ps_prop.setShuffle();
    ps_prop.setDeflate(compression);

    _ps_dataset = _file.createDataSet("/PhaseSpace/data",_ps_datatype,
                                       ps_dataspace,ps_prop);

    const double ps_AmperePerSigma2 = ps->current;
    const double ps_CoulombPerSigma2 = ps->charge;

    _ps_dataset.createAttribute("AmperePerNBLPerNES",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &ps_AmperePerSigma2);
    _ps_dataset.createAttribute("CoulombPerNBLPerNES",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &ps_CoulombPerSigma2);
    }


    if (imp != nullptr) {
        // save Impedance
        _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_f", "/Impedance/axis0" );


        if (std::is_same<vfps::impedance_t,std::complex<float>>::value) {
            imp_datatype = H5::PredType::IEEE_F32LE;
        } else if (std::is_same<vfps::impedance_t,std::complex<fixp64>>::value) {
            imp_datatype = H5::PredType::STD_I64LE;
        } else if (std::is_same<vfps::impedance_t,std::complex<double>>::value) {
            imp_datatype = H5::PredType::IEEE_F64LE;
        }

        H5::DataSpace imp_dataspace(imp_rank,&imp_size,&imp_size);

        const hsize_t imp_chunkdims = std::min(hsize_t(4096),imp_size);

        H5::DSetCreatPropList imp_prop;
        imp_prop.setChunk(imp_rank,&imp_chunkdims);
        imp_prop.setShuffle();
        imp_prop.setDeflate(compression);

        _file.openGroup("Impedance/data").createAttribute
            ("Ohm",H5::PredType::IEEE_F64LE,
             H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&(imp->factor4Ohms));
        imp_dataset_real = _file.createDataSet("/Impedance/data/real",imp_datatype,
                                                imp_dataspace,imp_prop);
        imp_dataset_imag = _file.createDataSet("/Impedance/data/imag",imp_datatype,
                                                imp_dataspace,imp_prop);


        std::vector<csrpower_t> imp_real;
        std::vector<csrpower_t> imp_imag;
        imp_real.reserve(imp_size);
        imp_imag.reserve(imp_size);
        for (const impedance_t& z : imp->impedance()) {
            imp_real.push_back(z.real());
            imp_imag.push_back(z.imag());
        }
        imp_dataset_real.write(imp_real.data(),imp_datatype);
        imp_dataset_imag.write(imp_imag.data(),imp_datatype);
    }

    // get ready to save SourceMap
    {
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/SourceMap/axis0" );
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/SourceMap/axis1" );
    _file.link(H5L_TYPE_SOFT, "/Info/AxisValues_E", "/SourceMap/axis2" );

    _sm_datatype = H5::PredType::IEEE_F32LE;

    const std::array<hsize_t,_sm_rank> sm_maxdims
            = {{H5S_UNLIMITED,_nBunches,_psSizeX,_psSizeY}};

    H5::DataSpace sm_dataspace(_sm_rank,_sm_dims.data(),sm_maxdims.data());


    const std::array<hsize_t,_sm_rank> sm_chunkdims
        = {{64U,1U,std::min(64U,_psSizeX),std::min(64U,_psSizeY)}};

    H5::DSetCreatPropList sm_prop;
    sm_prop.setChunk(_sm_rank,sm_chunkdims.data());
    sm_prop.setShuffle();
    sm_prop.setDeflate(compression);

    _sm_dataset_x = _file.createDataSet("/SourceMap/data/x",_sm_datatype,
                                       sm_dataspace,sm_prop);
    _sm_dataset_y = _file.createDataSet("/SourceMap/data/y",_sm_datatype,
                                       sm_dataspace,sm_prop);
    }

    if (wfm != nullptr ) {
        // save Wake Function

        if (std::is_same<vfps::meshaxis_t,float>::value) {
            wf_datatype = H5::PredType::IEEE_F32LE;
        } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            wf_datatype = H5::PredType::STD_I64LE;
        } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            wf_datatype = H5::PredType::IEEE_F64LE;
        }

        H5::DataSpace wf_dataspace(wf_rank,&wf_size,&wf_size);

        const hsize_t wf_chunkdims = std::min(hsize_t(4096),wf_size);

        H5::DSetCreatPropList wf_prop;
        wf_prop.setChunk(wf_rank,&wf_chunkdims);
        wf_prop.setShuffle();
        wf_prop.setDeflate(compression);

        wf_dataset = _file.createDataSet("/WakeFunction/data",wf_datatype,
                                          wf_dataspace,wf_prop);

        wf_dataset.write(wfm->getWakeFunction(),wf_datatype);
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
    H5::DataSpace* filespace;
    H5::DataSpace* memspace;


    if (fullspectrum) {
        // append CSR Field
        std::array<hsize_t,csr_rank> csr_offset
                        = {{csr_dims[0],0}};
        const std::array<hsize_t,csr_rank> csr_ext
                        = {{1,csr_dims[1]}};
        csr_dims[0]++;
        csr_dataset.extend(csr_dims.data());
        filespace = new H5::DataSpace(csr_dataset.getSpace());
        filespace->selectHyperslab(H5S_SELECT_SET, csr_ext.data(), csr_offset.data());
        memspace = new H5::DataSpace(csr_rank,csr_ext.data(),nullptr);
        csr_dataset.write(ef->getCSRSpectrum(), csr_datatype,*memspace, *filespace);
        delete memspace;
        delete filespace;
    }

    // append CSR Power
    std::array<hsize_t,csri_rank> csri_offset({{csri_dims[0]}});
    const std::array<hsize_t,csri_rank> csri_ext
                         = {{1,csri_dims[1]}};
    csri_dims[0]++;
    csri_dataset.extend(csri_dims.data());
    filespace = new H5::DataSpace(csri_dataset.getSpace());
    filespace->selectHyperslab(H5S_SELECT_SET, csri_ext.data(), csri_offset.data());
    memspace = new H5::DataSpace(csri_rank,csri_ext.data(),nullptr);
    csrpower_t csrpower = ef->getCSRPower();
    csri_dataset.write(&csrpower, csri_datatype,*memspace, *filespace);
    delete memspace;
    delete filespace;
}

void vfps::HDF5File::appendRFKicks(
        const std::vector<std::array<vfps::meshaxis_t,2>> kicks)
{
    std::array<hsize_t,drfk_rank> drfk_offset = {{drfk_dims[0],0}};
    const std::array<hsize_t,drfk_rank> drfk_ext = {{kicks.size(),2}};
    drfk_dims[0] += kicks.size();
    drfk_dataset.extend(drfk_dims.data());
    H5::DataSpace filespace(drfk_dataset.getSpace());
    filespace.selectHyperslab( H5S_SELECT_SET, drfk_ext.data()
                              , drfk_offset.data());
    H5::DataSpace memspace( drfk_rank, drfk_ext.data(), nullptr);
    drfk_dataset.write(kicks.data(), drfk_datatype,memspace, filespace);
}

void vfps::HDF5File::appendTracks(const PhaseSpace::Position *particles)
{
    std::array<hsize_t,pt_rank> pt_offset
            = {{pt_dims[0],0,0}};
    const std::array<hsize_t,pt_rank> pt_ext
            = {{1,pt_particles,2}};
    pt_dims[0]++;
    pt_dataset.extend(pt_dims.data());
    H5::DataSpace* filespace = new H5::DataSpace(pt_dataset.getSpace());
    filespace->selectHyperslab(H5S_SELECT_SET, pt_ext.data(), pt_offset.data());
    H5::DataSpace* memspace = new H5::DataSpace(pt_rank,pt_ext.data(),nullptr);
    pt_dataset.write(particles, pt_datatype,*memspace, *filespace);
    delete memspace;
    delete filespace;
}

void vfps::HDF5File::appendSourceMap(const PhaseSpace::Position *allpos)
{
    H5::DataSpace* filespace = nullptr;
    H5::DataSpace* memspace = nullptr;

    float* data_x = new float[_psSizeX*_psSizeX];
    float* data_y = new float[_psSizeX*_psSizeX];
    // get deltas from positions
    for (meshindex_t x=0; x<_psSizeX; x++) {
        for (meshindex_t y=0; y<_psSizeX; y++) {
            data_x[x*_psSizeX+y] = allpos[x*_psSizeX+y].x-x;
            data_y[x*_psSizeX+y] = allpos[x*_psSizeX+y].y-y;
        }
    }

   // append SourceMap
   std::array<hsize_t,_sm_rank> sm_offset
           = {{_sm_dims[0],0,0}};
   const std::array<hsize_t,_sm_rank> sm_ext
           = {{1,_psSizeX,_psSizeY}};
   _sm_dims[0]++;
   memspace = new H5::DataSpace(_sm_rank,sm_ext.data(),nullptr);

   _sm_dataset_x.extend(_sm_dims.data());
   filespace = new H5::DataSpace(_sm_dataset_x.getSpace());
   filespace->selectHyperslab(H5S_SELECT_SET, sm_ext.data(), sm_offset.data());
   _sm_dataset_x.write(data_x, _sm_datatype, *memspace, *filespace);

   _sm_dataset_y.extend(_sm_dims.data());
   filespace = new H5::DataSpace(_sm_dataset_y.getSpace());
   filespace->selectHyperslab(H5S_SELECT_SET, sm_ext.data(), sm_offset.data());
   _sm_dataset_y.write(data_y, _sm_datatype, *memspace, *filespace);

   delete memspace;
   delete filespace;
   delete [] data_x;
   delete [] data_y;
}

void vfps::HDF5File::append(const PhaseSpace& ps,
                            const timeaxis_t t,
                            const AppendType at)
{
    H5::DataSpace* filespace = nullptr;
    H5::DataSpace* memspace = nullptr;

    if ( at == AppendType::All ||
         at == AppendType::PhaseSpace) {
        // append to TimeAxis
        hsize_t ps_ta_offset = _ps_ta_dims;
        const hsize_t ps_ta_ext = 1;
        _ps_ta_dims++;
        _ps_ta_dataset.extend(&_ps_ta_dims);
        filespace = new H5::DataSpace(_ps_ta_dataset.getSpace());
        filespace->selectHyperslab(H5S_SELECT_SET, &ps_ta_ext, &ps_ta_offset);
        memspace = new H5::DataSpace(timeaxis.rank,&ps_ta_ext,nullptr);
        _ps_ta_dataset.write(&t, timeaxis.datatype,*memspace, *filespace);
        delete memspace;
        delete filespace;

        // append PhaseSpace
        std::array<hsize_t,_ps_rank> ps_offset
                = {{_ps_dims[0],0,0}};
        const std::array<hsize_t,_ps_rank> ps_ext
                = {{1,_psSizeX,_psSizeX}};
        _ps_dims[0]++;
        _ps_dataset.extend(_ps_dims.data());
        filespace = new H5::DataSpace(_ps_dataset.getSpace());
        memspace = new H5::DataSpace(_ps_rank,ps_ext.data(),nullptr);
        filespace->selectHyperslab(H5S_SELECT_SET, ps_ext.data(), ps_offset.data());
        _ps_dataset.write(ps.getData(), _ps_datatype, *memspace, *filespace);
        delete memspace;
        delete filespace;
    }

    if (at != AppendType::PhaseSpace) {
        _appendData(timeaxis,&t);
        _appendData(bunchprofile,ps.getProjection(0));
        _appendData(bunchlength,ps.getBunchLength());
        {
        auto mean_q = ps.getMoment(0,0);
        _appendData(bunchposition,&mean_q);
        }

        // append energy profile
        {
        std::array<hsize_t,ep_rank> ep_offset
                        = {{ep_dims[0],0}};
        const std::array<hsize_t,ep_rank> ep_ext
                        = {{1,_psSizeX}};
        ep_dims[0]++;
        ep_dataset.extend(ep_dims.data());
        H5::DataSpace filespace(ep_dataset.getSpace());
        filespace.selectHyperslab(H5S_SELECT_SET, ep_ext.data(), ep_offset.data());
        H5::DataSpace memspace(ep_rank,ep_ext.data(),nullptr);
        ep_dataset.write(ps.getProjection(1), ep_datatype,memspace,filespace);
        }

        // append energy spread
        {
        std::array<hsize_t,es_rank> es_offset({{es_dims[0]}});
        const std::array<hsize_t,es_rank> es_ext{{ 1, es_dims[1] }};
        es_dims[0]++;
        es_dataset.extend(es_dims.data());
        H5::DataSpace filespace(es_dataset.getSpace());
        filespace.selectHyperslab(H5S_SELECT_SET, es_ext.data(), es_offset.data());
        H5::DataSpace memspace(es_rank,es_ext.data(),nullptr);
        auto energyspread = ps.getEnergySpread();
        es_dataset.write(&energyspread, es_datatype,memspace, filespace);
        }

        // append BunchPopulation
        {
        double bc = ps.getIntegral();
        _appendData(bunchpopulation,&bc);
        }
    }
}

void vfps::HDF5File::append(const WakeKickMap* wkm)
{
    // append WakePotential
    std::array<hsize_t,wp_rank> wp_offset
            = {{wp_dims[0],0}};
    const std::array<hsize_t,wp_rank> wp_ext
            = {{1,_psSizeX}};
    wp_dims[0]++;
    wp_dataset.extend(wp_dims.data());
    H5::DataSpace* filespace = new H5::DataSpace(wp_dataset.getSpace());
    filespace->selectHyperslab(H5S_SELECT_SET, wp_ext.data(), wp_offset.data());
    H5::DataSpace* memspace = new H5::DataSpace(wp_rank,wp_ext.data(),nullptr);
    wp_dataset.write(wkm->getForce(), wp_datatype,*memspace, *filespace);
    delete memspace;
    delete filespace;
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

    hsize_t ps_dims[3];
    ps_space.getSimpleExtentDims( ps_dims, nullptr );

    meshindex_t ps_size = ps_dims[1];
    use_step = (ps_dims[0]+use_step)%ps_dims[0];
    const std::array<hsize_t,_ps_rank> ps_offset =
        {{static_cast<hsize_t>(use_step),0,0}};
    const std::array<hsize_t,_ps_rank> ps_ext = {{1,ps_size,ps_size}};
    H5::DataSpace memspace(_ps_rank,ps_ext.data(),nullptr);
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
void vfps::HDF5File::_appendData( DatasetInfo<rank>& ds
                                , const datatype* const data)
{
    std::array<hsize_t,ds.rank> offset{{ds.dims[0]}};
    auto ext = ds.dims;
    ext[0] = 1;
    ds.dims[0]++;
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
