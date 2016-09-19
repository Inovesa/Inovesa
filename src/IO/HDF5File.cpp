/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2016: Karlsruhe Institute of Technology                 *
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

vfps::HDF5File::HDF5File(const std::string filename,
                         const PhaseSpace* ps,
                         const ElectricField* ef,
                         const Impedance* imp,
                         const WakeFunctionMap* wfm,
                         const size_t nparticles,
                         const double t_sync) :
    _file( nullptr ),
    fname( filename ),
    ta_dims( 0 ),
    bc_dims( 0 ),
    bp_dims( {{ 0, ps->nMeshCells(0) }} ),
    bl_dims( 0 ),
    qb_dims( 0 ),
    ep_dims( {{ 0, ps->nMeshCells(0) }} ),
    es_dims( 0 ),
    pt_dims( {{ 0, nparticles, 2 }} ),
    pt_particles( nparticles ),
    wp_dims( {{ 0, ps->nMeshCells(0) }} ),
    csr_dims( {{ 0, ef->getNMax()/static_cast<size_t>(2) }} ),
    maxn( ef->getNMax()/static_cast<size_t>(2) ),
    csri_dims( 0 ),
    ps_dims( {{ 0, ps->nMeshCells(0), ps->nMeshCells(1) }} ),
    ps_size( ps->nMeshCells(0) ),
    imp_size( imp->nFreqs()/2 ),
    wf_size( 2*ps_size )
{
    _file = new H5::H5File(filename,H5F_ACC_TRUNC);

    _file->createGroup("Info");
    // save Values of Phase Space Axis
    if (std::is_same<vfps::meshaxis_t,float>::value) {
            axps_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            axps_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            axps_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,axps_rank> psa_dims
                    = {{ ps_size }};
    const std::array<hsize_t,axps_rank> psa_maxdims
                    = {{ ps_size }};

    H5::DataSpace ax0ps_dataspace(axps_rank,psa_dims.data(),
                                        psa_maxdims.data());
    H5::DataSpace ax1ps_dataspace(axps_rank,psa_dims.data(),
                                        psa_maxdims.data());


    const hsize_t psa_chunkdims = std::min(2048U,ps_size);
    axps_prop.setChunk(axps_rank,&psa_chunkdims);
    axps_prop.setShuffle();
    axps_prop.setDeflate(compression);

    ax0ps_dataset = _file->createDataSet("/Info/AxisValues_z",axps_datatype,
                                         ax0ps_dataspace,axps_prop);
    const double ax0scale0 = ps->getScale(0);
    ax0ps_dataset.createAttribute("Factor4Meters",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax0scale0);
    const double ax0scale1 = ax0scale0/physcons::c;
    ax0ps_dataset.createAttribute("Factor4Seconds",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax0scale1);

    ax1ps_dataset = _file->createDataSet("/Info/AxisValues_E",axps_datatype,
                                         ax1ps_dataspace,axps_prop);
    const double ax1scale = ps->getScale(1);
    ax1ps_dataset.createAttribute("Factor4ElectronVolts"
                                  ,H5::PredType::IEEE_F64LE,
                                  H5::DataSpace()).write(
          H5::PredType::IEEE_F64LE,&ax1scale);
    ax0ps_dataset.write(ps->getAxis(0)->data(),axps_datatype);
    ax1ps_dataset.write(ps->getAxis(0)->data(),axps_datatype);


    // create group for parameters
    _file->createGroup("/Info/Parameters");

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
    axfreq_prop.setChunk(axps_rank,&axfreq_chunkdims);
    axfreq_prop.setShuffle();
    axfreq_prop.setDeflate(compression);

    axfreq_dataset = _file->createDataSet("/Info/AxisValues_f",axfreq_datatype,
                                          axfreq_dataspace,axfreq_prop);
    const double axfreqscale = imp->getRuler()->scale();
    axfreq_dataset.createAttribute("Factor4Hertz",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&axfreqscale);
    axfreq_dataset.write(imp->getRuler()->data(),axfreq_datatype);


    // get ready to save TimeAxis
    ta_datatype = H5::PredType::IEEE_F64LE;

    hsize_t ta_maxdims = H5S_UNLIMITED;

    H5::DataSpace ta_dataspace(ta_rank,&ta_dims,&ta_maxdims);


    const hsize_t ta_chunkdims = 256;
    ta_prop.setChunk(ta_rank,&ta_chunkdims);
    ta_prop.setShuffle();
    ta_prop.setDeflate(compression);

    ta_dataset = _file->createDataSet("/Info/AxisValues_t",ta_datatype,
                                      ta_dataspace,ta_prop);
    const double axtimescale = t_sync;
    ta_dataset.createAttribute("Factor4Seconds",H5::PredType::IEEE_F64LE,
                H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&axtimescale);

    // get ready to save BunchCharge
    _file->createGroup("BunchPopulation");
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_t","/BunchPopulation/axis0");
    bc_datatype = H5::PredType::IEEE_F64LE;

    hsize_t bc_maxdims = H5S_UNLIMITED;

     H5::DataSpace bc_dataspace(bc_rank,&bc_dims,&bc_maxdims);

    const hsize_t bc_chunkdims = 256;
    bc_prop.setChunk(bc_rank,&bc_chunkdims);
    bc_prop.setShuffle();
    bc_prop.setDeflate(compression);

    bc_dataset = _file->createDataSet("/BunchPopulation/data",bc_datatype,
                                      bc_dataspace,bc_prop);
    bc_dataset.createAttribute("Factor4Ampere",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&(ps->current));
    bc_dataset.createAttribute("Factor4Coulomb",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&(ps->charge));

    // get ready to save BunchProfiles
    _file->createGroup("BunchProfile");
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/BunchProfile/axis0" );
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/BunchProfile/axis1" );

    if (std::is_same<vfps::integral_t,float>::value) {
            bp_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::integral_t,fixp64>::value) {
            bp_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::integral_t,double>::value) {
            bp_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,bp_rank> bp_maxdims
                    = {{H5S_UNLIMITED,ps_size}};

    H5::DataSpace bp_dataspace(bp_rank,bp_dims.data(),bp_maxdims.data());

    const std::array<hsize_t,bp_rank> bp_chunkdims
        = {{64U,std::min(2048U,ps_size)}};
    bp_prop.setChunk(bp_rank,bp_chunkdims.data());
    bp_prop.setShuffle();
    bp_prop.setDeflate(compression);

    const double bp_factor4ampere = ps->getAxis(0)->delta()*ps->current;
    const double bp_factor4coulomb = ps->getAxis(0)->delta()*ps->charge;

    bp_dataset = _file->createDataSet("/BunchProfile/data",bp_datatype,
                                      bp_dataspace,bp_prop);
    bp_dataset.createAttribute("Factor4Ampere",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&bp_factor4ampere);
    bp_dataset.createAttribute("Factor4Coulomb",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&bp_factor4coulomb);

    // get ready to save BunchLength
    _file->createGroup("BunchLength");
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/BunchLength/axis0" );
    if (std::is_same<vfps::meshaxis_t,float>::value) {
            bl_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            bl_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            bl_datatype = H5::PredType::IEEE_F64LE;
    }

    hsize_t bl_maxdims = H5S_UNLIMITED;

    H5::DataSpace bl_dataspace(bl_rank,&bl_dims,&bl_maxdims);

    const hsize_t bl_chunkdims = 256;
    bl_prop.setChunk(bl_rank,&bl_chunkdims);
    bl_prop.setShuffle();
    bl_prop.setDeflate(compression);

    bl_dataset = _file->createDataSet("/BunchLength/data",bl_datatype,
                                      bl_dataspace,bl_prop);
    bl_dataset.createAttribute("Factor4Meters",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax0scale0);
    bl_dataset.createAttribute("Factor4Seconds",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax0scale1);

    // get ready to save BunchPosition
    _file->createGroup("BunchPosition");
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/BunchPosition/axis0" );
    if (std::is_same<vfps::meshaxis_t,float>::value) {
            qb_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            qb_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            qb_datatype = H5::PredType::IEEE_F64LE;
    }

    hsize_t qb_maxdims = H5S_UNLIMITED;

    H5::DataSpace qb_dataspace(qb_rank,&qb_dims,&qb_maxdims);

    const hsize_t qb_chunkdims = 256;
    qb_prop.setChunk(qb_rank,&qb_chunkdims);
    qb_prop.setShuffle();
    qb_prop.setDeflate(compression);

    qb_dataset = _file->createDataSet("/BunchPosition/data",qb_datatype,
                                      qb_dataspace,qb_prop);
    qb_dataset.createAttribute("Factor4Meters",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax0scale0);
    qb_dataset.createAttribute("Factor4Seconds",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax0scale1);



    // get ready to save EnergyProfiles
    _file->createGroup("EnergyProfile");
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_E", "/EnergyProfile/axis0" );
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/EnergyProfile/axis1" );

    if (std::is_same<vfps::integral_t,float>::value) {
            ep_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::integral_t,fixp64>::value) {
            ep_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::integral_t,double>::value) {
            ep_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,ep_rank> ep_maxdims
                    = {{H5S_UNLIMITED,ps_size}};

    H5::DataSpace ep_dataspace(ep_rank,ep_dims.data(),ep_maxdims.data());

    const std::array<hsize_t,ep_rank> ep_chunkdims
        = {{64U,std::min(2048U,ps_size)}};
    ep_prop.setChunk(ep_rank,ep_chunkdims.data());
    ep_prop.setShuffle();
    ep_prop.setDeflate(compression);


    const double ep_factor4ampere = ps->getAxis(1)->delta()*ps->current;
    const double ep_factor4coulomb = ps->getAxis(1)->delta()*ps->charge;

    ep_dataset = _file->createDataSet("/EnergyProfile/data",ep_datatype,
                                      ep_dataspace,ep_prop);

    ep_dataset.createAttribute("Factor4Ampere",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ep_factor4ampere);
    ep_dataset.createAttribute("Factor4Coulomb",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ep_factor4coulomb);

    // get ready to save Energy Spread
    _file->createGroup("EnergySpread");
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/EnergySpread/axis0" );
    if (std::is_same<vfps::meshaxis_t,float>::value) {
            es_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            es_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            es_datatype = H5::PredType::IEEE_F64LE;
    }

    hsize_t es_maxdims = H5S_UNLIMITED;

    H5::DataSpace es_dataspace(es_rank,&es_dims,&es_maxdims);

    const hsize_t es_chunkdims = 256;
    es_prop.setChunk(es_rank,&es_chunkdims);
    es_prop.setShuffle();
    es_prop.setDeflate(compression);

    es_dataset = _file->createDataSet("/EnergySpread/data",es_datatype,
                                      es_dataspace,es_prop);
    es_dataset.createAttribute("Factor4ElectronVolts",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ax1scale);

    // get ready to save particles from (pseudo-) tracking
    _file->createGroup("Particles");
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/Particles/axis0" );
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
    pt_prop.setChunk(pt_rank,pt_chunkdims.data());
    pt_prop.setShuffle();
    pt_prop.setDeflate(compression);

    pt_dataset = _file->createDataSet("/Particles/data",pt_datatype,
                                      pt_dataspace,pt_prop);

    // get ready to save WakePotential
    _file->createGroup("WakePotential");
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/WakePotential/axis0" );

    if (std::is_same<vfps::meshaxis_t,float>::value) {
        wp_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
        wp_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::meshaxis_t,double>::value) {
        wp_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,wp_rank> wp_maxdims
            = {{H5S_UNLIMITED,ps_size}};

    H5::DataSpace wp_dataspace(wp_rank,wp_dims.data(),wp_maxdims.data());

    const std::array<hsize_t,wp_rank> wp_chunkdims
            = {{64U,std::min(2048U,ps_size)}};
    wp_prop.setChunk(wp_rank,wp_chunkdims.data());
    wp_prop.setShuffle();
    wp_prop.setDeflate(compression);

    wp_dataset = _file->createDataSet("/WakePotential/data",wp_datatype,
                                      wp_dataspace,wp_prop);
    wp_dataset.createAttribute("Factor4Volts",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,
                                   &(ef->volts));

    // get ready to save CSR Data
    _file->createGroup("CSR/");
    // get ready to save CSR Spectrum
    _file->createGroup("CSR/Spectrum");
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/CSR/Spectrum/axis0" );
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_f", "/CSR/Spectrum/axis1" );

    if (std::is_same<vfps::csrpower_t,float>::value) {
        csr_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::csrpower_t,double>::value) {
        csr_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,csr_rank> csr_maxdims
        = {{H5S_UNLIMITED,csr_dims[1]}};

    H5::DataSpace csr_dataspace(csr_rank,csr_dims.data(),csr_maxdims.data());


    const std::array<hsize_t,csr_rank> csr_chunkdims
        = {{64U,std::min(static_cast<size_t>(2048),maxn)}};
    csr_prop.setChunk(csr_rank,csr_chunkdims.data());
    csr_prop.setShuffle();
    csr_prop.setDeflate(compression);

    const double csr_factor4watts = std::pow(bp_factor4ampere,2)
                                  * imp->factor4Ohms;

    csr_dataset = _file->createDataSet("/CSR/Spectrum/data",csr_datatype,
                                       csr_dataspace,csr_prop);
    csr_dataset.createAttribute("Factor4Watts",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE, &csr_factor4watts);

    // get ready to save CSR Intensity
    _file->createGroup("CSR/Intensity");
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_t", "/CSR/Intensity/axis0" );

    if (std::is_same<vfps::csrpower_t,float>::value) {
        csri_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::csrpower_t,double>::value) {
        csri_datatype = H5::PredType::IEEE_F64LE;
    }

    const hsize_t csrp_maxdims = H5S_UNLIMITED;

    H5::DataSpace csrp_dataspace(csri_rank,&csri_dims,&csrp_maxdims);


    const hsize_t csrp_chunkdims = std::min(2048U,ps_size);
    csri_prop.setChunk(csri_rank,&csrp_chunkdims);
    csri_prop.setShuffle();
    csri_prop.setDeflate(compression);

    csri_dataset = _file->createDataSet("/CSR/Intensity/data",csri_datatype,
                                        csrp_dataspace,csri_prop);

    csri_dataset.createAttribute("Factor4Watts",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE, &csr_factor4watts);

    // get ready to save PhaseSpace
    _file->createGroup("PhaseSpace");
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_z", "/PhaseSpace/axis0" );
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_E", "/PhaseSpace/axis1" );

    if (std::is_same<vfps::meshdata_t,float>::value) {
        ps_datatype = H5::PredType::IEEE_F32LE;
    #if FXP_FRACPART < 31
    } else if (std::is_same<vfps::meshdata_t,fixp32>::value) {
        ps_datatype = H5::PredType::STD_I32LE;
    #endif
    } else if (std::is_same<vfps::meshdata_t,fixp64>::value) {
        ps_datatype = H5::PredType::STD_I64LE;
    }else if (std::is_same<vfps::meshdata_t,double>::value) {
        ps_datatype = H5::PredType::IEEE_F64LE;
    }

    const std::array<hsize_t,ps_rank> ps_maxdims
            = {{H5S_UNLIMITED,ps_size,ps_size}};

    H5::DataSpace ps_dataspace(ps_rank,ps_dims.data(),ps_maxdims.data());


    const std::array<hsize_t,ps_rank> ps_chunkdims
        = {{64U,std::min(64U,ps_size),std::min(64U,ps_size)}};
    ps_prop.setChunk(ps_rank,ps_chunkdims.data());
    ps_prop.setShuffle();
    ps_prop.setDeflate(compression);

    _ps_dataset = _file->createDataSet("/PhaseSpace/data",ps_datatype,
                                       ps_dataspace,ps_prop);

    const double ps_factor4ampere = ps->getAxis(0)->delta()
                                  * ps->getAxis(1)->delta()
                                  * ps->current;
    const double ps_factor4coulomb = ps->getAxis(0)->delta()
                                   * ps->getAxis(1)->delta()
                                   * ps->charge;

    _ps_dataset.createAttribute("Factor4Ampere",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ps_factor4ampere);
    _ps_dataset.createAttribute("Factor4Coulomb",H5::PredType::IEEE_F64LE,
            H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&ps_factor4coulomb);

    // save Impedance
    _file->createGroup("Impedance");
    _file->link(H5L_TYPE_SOFT, "/Info/AxisValues_f", "/Impedance/axis0" );

    if (std::is_same<vfps::impedance_t,std::complex<float>>::value) {
        imp_datatype = H5::PredType::IEEE_F32LE;
    } else if (std::is_same<vfps::impedance_t,std::complex<fixp64>>::value) {
        imp_datatype = H5::PredType::STD_I64LE;
    } else if (std::is_same<vfps::impedance_t,std::complex<double>>::value) {
        imp_datatype = H5::PredType::IEEE_F64LE;
    }

    H5::DataSpace imp_dataspace(imp_rank,&imp_size,&imp_size);

    const hsize_t imp_chunkdims = std::min(hsize_t(4096),imp_size);
    imp_prop.setChunk(imp_rank,&imp_chunkdims);
    imp_prop.setShuffle();
    imp_prop.setDeflate(compression);

    _file->createGroup("Impedance/data").createAttribute
        ("Factor4Ohms",H5::PredType::IEEE_F64LE,
         H5::DataSpace()).write(H5::PredType::IEEE_F64LE,&(imp->factor4Ohms));
    imp_dataset_real = _file->createDataSet("/Impedance/data/real",imp_datatype,
                                            imp_dataspace,imp_prop);
    imp_dataset_imag = _file->createDataSet("/Impedance/data/imag",imp_datatype,
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

    if (wfm != nullptr ) {
        // save Wake Function
        _file->createGroup("WakeFunction");

        if (std::is_same<vfps::meshaxis_t,float>::value) {
            wf_datatype = H5::PredType::IEEE_F32LE;
        } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
            wf_datatype = H5::PredType::STD_I64LE;
        } else if (std::is_same<vfps::meshaxis_t,double>::value) {
            wf_datatype = H5::PredType::IEEE_F64LE;
        }

        H5::DataSpace wf_dataspace(wf_rank,&wf_size,&wf_size);

        const hsize_t wf_chunkdims = std::min(hsize_t(4096),wf_size);
        wf_prop.setChunk(wf_rank,&wf_chunkdims);
        wf_prop.setShuffle();
        wf_prop.setDeflate(compression);

        wf_dataset = _file->createDataSet("/WakeFunction/data",wf_datatype,
                                          wf_dataspace,wf_prop);

        wf_dataset.write(wfm->getWakeFunction(),wf_datatype);
    }

    // save Inovesa version
    std::array<hsize_t,1> version_dims {{3}};
    H5::DataSpace version_dspace(1,version_dims.data(),version_dims.data());
    H5::DataSet version_dset = _file->createDataSet
                    ("/Info/INOVESA_v", H5::PredType::STD_I32LE,version_dspace);
    std::array<int32_t,3> version {{INOVESA_VERSION_RELEASE,
                                    INOVESA_VERSION_MINOR,
                                    INOVESA_VERSION_FIX}};
    version_dset.write(version.data(),H5::PredType::NATIVE_INT);
}

vfps::HDF5File::~HDF5File()
{
    delete _file;
}

void vfps::HDF5File::addParameterToGroup(std::string groupname,
                                         std::string paramname,
                                         H5::PredType type,
                                         void* data)
{
    H5::Group group = _file->openGroup(groupname);
    group.createAttribute(paramname,type, H5::DataSpace()).write(type,data);
}

void vfps::HDF5File::append(const ElectricField* ef, bool fullspectrum)
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

    // append CSR Intensity
    hsize_t csrp_offset = csri_dims;
    const hsize_t csrp_ext = 1;
    csri_dims++;
    csri_dataset.extend(&csri_dims);
    filespace = new H5::DataSpace(csri_dataset.getSpace());
    filespace->selectHyperslab(H5S_SELECT_SET, &csrp_ext, &csrp_offset);
    memspace = new H5::DataSpace(csri_rank,&csrp_ext,nullptr);
    csrpower_t csrpower = ef->getCSRPower();
    csri_dataset.write(&csrpower, csri_datatype,*memspace, *filespace);
    delete memspace;
    delete filespace;
}

void vfps::HDF5File::append(const PhaseSpace::Position *particles)
{
    // append WakePotential
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

void vfps::HDF5File::append(const PhaseSpace* ps, const AppendType at)
{
    H5::DataSpace* filespace;
    H5::DataSpace* memspace;

    if ( at == AppendType::All ||
         at == AppendType::PhaseSpace) {
        // append PhaseSpace
        std::array<hsize_t,ps_rank> ps_offset
                = {{ps_dims[0],0,0}};
        const std::array<hsize_t,ps_rank> ps_ext
                = {{1,ps_size,ps_size}};
        ps_dims[0]++;
        _ps_dataset.extend(ps_dims.data());
        filespace = new H5::DataSpace(_ps_dataset.getSpace());
        memspace = new H5::DataSpace(ps_rank,ps_ext.data(),nullptr);
        filespace->selectHyperslab(H5S_SELECT_SET, ps_ext.data(), ps_offset.data());
        _ps_dataset.write(ps->getData(), ps_datatype, *memspace, *filespace);
        delete memspace;
        delete filespace;
    }

    if (at != AppendType::PhaseSpace) {
        // append BunchProfile
        std::array<hsize_t,bp_rank> bp_offset
                        = {{bp_dims[0],0}};
        const std::array<hsize_t,bp_rank> bp_ext
                        = {{1,ps_size}};
        bp_dims[0]++;
        bp_dataset.extend(bp_dims.data());
        filespace = new H5::DataSpace(bp_dataset.getSpace());
        filespace->selectHyperslab(H5S_SELECT_SET, bp_ext.data(), bp_offset.data());
        memspace = new H5::DataSpace(bp_rank,bp_ext.data(),nullptr);
        bp_dataset.write(ps->getProjection(0), bp_datatype,*memspace, *filespace);
        delete memspace;
        delete filespace;

        // append Length
        hsize_t bl_offset = bl_dims;
        const hsize_t bl_ext = 1;
        bl_dims++;
        bl_dataset.extend(&bl_dims);
        filespace = new H5::DataSpace(bl_dataset.getSpace());
        filespace->selectHyperslab(H5S_SELECT_SET, &bl_ext, &bl_offset);
        memspace = new H5::DataSpace(bl_rank,&bl_ext,nullptr);
        meshaxis_t bunchlength = ps->getMoment(0,1);
        bl_dataset.write(&bunchlength, bl_datatype,*memspace, *filespace);
        delete memspace;
        delete filespace;

        // append Position
        hsize_t qb_offset = qb_dims;
        const hsize_t qb_ext = 1;
        qb_dims++;
        qb_dataset.extend(&qb_dims);
        filespace = new H5::DataSpace(qb_dataset.getSpace());
        filespace->selectHyperslab(H5S_SELECT_SET, &qb_ext, &qb_offset);
        memspace = new H5::DataSpace(qb_rank,&qb_ext,nullptr);
        meshaxis_t bunchposition = ps->getMoment(0,0);
        qb_dataset.write(&bunchposition, qb_datatype,*memspace, *filespace);
        delete memspace;
        delete filespace;

        // append energy profile
        std::array<hsize_t,ep_rank> ep_offset
                        = {{ep_dims[0],0}};
        const std::array<hsize_t,ep_rank> ep_ext
                        = {{1,ps_size}};
        ep_dims[0]++;
        ep_dataset.extend(ep_dims.data());
        filespace = new H5::DataSpace(ep_dataset.getSpace());
        filespace->selectHyperslab(H5S_SELECT_SET, ep_ext.data(), ep_offset.data());
        memspace = new H5::DataSpace(ep_rank,ep_ext.data(),nullptr);
        ep_dataset.write(ps->getProjection(1), ep_datatype,*memspace, *filespace);
        delete memspace;
        delete filespace;

        // append energy spread
        hsize_t es_offset = es_dims;
        const hsize_t es_ext = 1;
        es_dims++;
        es_dataset.extend(&es_dims);
        filespace = new H5::DataSpace(es_dataset.getSpace());
        filespace->selectHyperslab(H5S_SELECT_SET, &es_ext, &es_offset);
        memspace = new H5::DataSpace(es_rank,&es_ext,nullptr);
        meshaxis_t energyspread = ps->getMoment(1,1);
        es_dataset.write(&energyspread, es_datatype,*memspace, *filespace);
        delete memspace;
        delete filespace;

        // append BunchPopulation
        hsize_t bc_offset = bc_dims;
        const hsize_t bc_ext = 1;
        bc_dims++;
        bc_dataset.extend(&bc_dims);
        filespace = new H5::DataSpace(bc_dataset.getSpace());
        filespace->selectHyperslab(H5S_SELECT_SET, &bc_ext, &bc_offset);
        memspace = new H5::DataSpace(bc_rank,&bc_ext,nullptr);
        double bc = ps->getIntegral();
        bc_dataset.write(&bc, bc_datatype,*memspace, *filespace);
        delete memspace;
        delete filespace;
    }
}

void vfps::HDF5File::append(const WakeKickMap* wkm)
{
    // append WakePotential
    std::array<hsize_t,wp_rank> wp_offset
            = {{wp_dims[0],0}};
    const std::array<hsize_t,wp_rank> wp_ext
            = {{1,ps_size}};
    wp_dims[0]++;
    wp_dataset.extend(wp_dims.data());
    H5::DataSpace* filespace = new H5::DataSpace(wp_dataset.getSpace());
    filespace->selectHyperslab(H5S_SELECT_SET, wp_ext.data(), wp_offset.data());
    H5::DataSpace* memspace = new H5::DataSpace(wp_rank,wp_ext.data(),nullptr);
    wp_dataset.write(wkm->getForce(), wp_datatype,*memspace, *filespace);
    delete memspace;
    delete filespace;
}

void vfps::HDF5File::appendTime(const double t)
{
    // append to TimeAxis
    hsize_t ta_offset = ta_dims;
    const hsize_t ta_ext = 1;
    ta_dims++;
    ta_dataset.extend(&ta_dims);
    H5::DataSpace* filespace;
    H5::DataSpace* memspace;
    filespace = new H5::DataSpace(ta_dataset.getSpace());
    filespace->selectHyperslab(H5S_SELECT_SET, &ta_ext, &ta_offset);
    memspace = new H5::DataSpace(ta_rank,&ta_ext,nullptr);
    ta_dataset.write(&t, ta_datatype,*memspace, *filespace);
    delete memspace;
    delete filespace;
}

vfps::PhaseSpace vfps::HDF5File::readPhaseSpace(std::string fname,
                                                meshaxis_t qmin,meshaxis_t qmax,
                                                meshaxis_t pmin,meshaxis_t pmax,
                                                double Qb, double Ib_unscaled,
                                                double bl, double dE)
{
    H5::DataType datatype;
    if (std::is_same<vfps::meshdata_t,float>::value) {
        datatype = H5::PredType::IEEE_F32LE;
    #if FXP_FRACPART < 31
    } else if (std::is_same<vfps::meshdata_t,fixp32>::value) {
        datatype = H5::PredType::STD_I32LE;
    #endif
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
    size_t ntimesteps = std::max(static_cast<hsize_t>(0),ps_dims[0]-1);
    const std::array<hsize_t,ps_rank> ps_offset = {{ntimesteps,0,0}};
    const std::array<hsize_t,ps_rank> ps_ext = {{1,ps_size,ps_size}};
    H5::DataSpace memspace(ps_rank,ps_ext.data(),nullptr);
    ps_space.selectHyperslab(H5S_SELECT_SET, ps_ext.data(), ps_offset.data());

    H5::DataType axistype;
    if (std::is_same<vfps::meshaxis_t,float>::value) {
        axistype = H5::PredType::IEEE_F32LE;
    #if FXP_FRACPART < 31
    } else if (std::is_same<vfps::meshaxis_t,fixp32>::value) {
        axistype = H5::PredType::STD_I32LE;
    #endif
    } else if (std::is_same<vfps::meshaxis_t,fixp64>::value) {
        axistype = H5::PredType::STD_I64LE;
    }else if (std::is_same<vfps::meshaxis_t,double>::value) {
        axistype = H5::PredType::IEEE_F64LE;
    }

    PhaseSpace ps(ps_size,qmin,qmax,pmin,pmax,Qb,Ib_unscaled,bl,dE);
    ps_dataset.read(ps.getData(), datatype, memspace, ps_space);

    return ps;
}

#endif // INOVESA_USE_HDF5
