// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#include "SM/WakeFunctionMap.hpp"

/**
 * @brief WakeFunctionMap
 * @param in
 * @param out
 * @param xsize
 * @param ysize
 * @param it
 * @param interpol_clamp
 * @param oclh
 *
 * @todo currently broken when used with CL/GL sharing
 */
vfps::WakeFunctionMap::WakeFunctionMap( std::shared_ptr<PhaseSpace> in
                                      , std::shared_ptr<PhaseSpace> out
                                      , const InterpolationType it
                                      , bool interpol_clamp
                                      , oclhptr_t oclh
                                      )
  : WakeKickMap( in,out,it,interpol_clamp,oclh
               #if INOVESA_USE_OPENCL ==1 and INOVESA_USE_OPENGL == 1
               , 0
               #endif // INOVESA_USE_OPENCL and INOVESA_USE_OPENGL
               )
  , _xaxis(Ruler<meshaxis_t>( 2*PhaseSpace::nx
                            , in->getMin(0)-in->length(0)/2
                            , in->getMax(0)+in->length(0)/2
                            , in->getScale(0)))
  , _wakefunction(new meshaxis_t[2*PhaseSpace::nx])
  , _wakesize(2*PhaseSpace::nx)
{
}

vfps::WakeFunctionMap::WakeFunctionMap( std::shared_ptr<PhaseSpace> in
                                      , std::shared_ptr<PhaseSpace> out
                                      , const std::string fname
                                      , const double sigmaE
                                      , const double E0
                                      , const double Ib
                                      , const double dt
                                      , const InterpolationType it
                                      , const bool interpol_clamp
                                      , oclhptr_t oclh
                                      )
  : WakeFunctionMap( in,out,it,interpol_clamp,oclh)
{
    /*  1e12*Ib*dt: wake file has charge in pC and for one revolution
     *  1/(ps->getDelta(1)*sigmaE*E0): eV -> pixels
     */
    _wakeFromFile(fname,1e12*Ib*dt/(in->getDelta(1)*E0*sigmaE));
}

vfps::WakeFunctionMap::WakeFunctionMap( std::shared_ptr<PhaseSpace> in
                                      , std::shared_ptr<PhaseSpace> out
                                      , const ElectricField* csr
                                      , const InterpolationType it
                                      , const bool interpol_clamp
                                      , oclhptr_t oclh
                                      )
  : WakeFunctionMap( in,out,it,interpol_clamp, oclh)
{
    std::copy_n(csr->getWakefunction(),2*PhaseSpace::nx,_wakefunction);
}

vfps::WakeFunctionMap::~WakeFunctionMap() noexcept
{
    delete [] _wakefunction;
    #if INOVESA_ENABLE_CLPROFILING == 1
    saveTimings("WakeFunctionMap");
    #endif // INOVESA_ENABLE_CLPROFILING
}

void vfps::WakeFunctionMap::update()
{
    #if INOVESA_USE_OPENCL
    if (_oclh) {
        _in->syncCLMem(OCLH::clCopyDirection::dev2cpu);
    }
    #endif // INOVESA_USE_OPENCL
    _in->integrate();
    integral_t charge = _in->getIntegral();
    auto density = _in->getProjection(0);
    for(uint32_t b=0; b<PhaseSpace::nb; b++) {
        for (unsigned int i=0;i<_xsize;i++) {
            _offset[b*_xsize+i] = 0;
            for (unsigned int j=0;j<_xsize;j++) {
                _offset[b*_xsize+i] += meshaxis_t(density[b][j]/charge
                                            *_wakefunction[_xsize+i-j]);
            }
        }
    }
    updateSM();
}

void vfps::WakeFunctionMap::_wakeFromFile(const std::string fname,
                                          const double scaling)
{
    std::ifstream ifs;
    ifs.open(fname);

    std::ofstream debug("debugwake.txt");

    std::vector<std::pair<meshaxis_t,double>> wake;

    while (ifs.good()) {
        double z,f;
        ifs >> z >> f;
        wake.push_back(std::pair<meshaxis_t,double>(z,f));
    }
    ifs.close();

    size_t x_other=1;
    size_t x_mine=0;
    std::array<interpol_t,2> qic;
    bool smallwake=false;

    if (_xaxis[0]*_xaxis.scale("Meter") < wake[0].first) {
        smallwake = true;

        // insert zeros left of the area in the wake file
        wake.insert(wake.begin(),std::pair<meshaxis_t,double>(wake[0].first-std::numeric_limits<meshaxis_t>::epsilon(),0));
        wake.insert(wake.begin(),std::pair<meshaxis_t,double>(_xaxis.min()*_xaxis.scale("Meter"),0));
    }

    // position at x_other should always be one left of the one of x_mine
    while (x_mine < _wakesize) {
        while (_xaxis[x_mine]*_xaxis.scale("Meter") > wake[x_other].first) {
            if (x_other+1 < wake.size()) {
                x_other++;
            } else {
                smallwake = true;
                // insert zeros right of the area in the wake file
                wake.push_back(std::pair<meshaxis_t,double>(wake.back().first+std::numeric_limits<meshaxis_t>::epsilon(),0));
                wake.push_back(std::pair<meshaxis_t,double>(_xaxis.max()*_xaxis.scale("Meter"),0));
            }
        }
        calcCoefficiants(qic.data(),
                         (wake[x_other].first-wake[x_other-1].first)/_xaxis.scale("Meter"),2);
        _wakefunction[x_mine]   = ( qic[0]*wake[x_other-1].second
                                  + qic[1]*wake[x_other  ].second )
                                * scaling;
        debug << x_mine << '\t'
              << _xaxis[x_mine]*_xaxis.scale("Meter") << '\t'
              << _wakefunction[x_mine] << '\t'
              << x_other << '\t'
              << wake[x_other-1].first << '\t'
              << wake[x_other  ].first << '\t'
              << wake[x_other  ].second << std::endl;
        x_mine++;
    }
    if (smallwake) {
        Display::printText("Warning: Given wake to small.");
    }
}
