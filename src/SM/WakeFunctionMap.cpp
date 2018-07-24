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

#include "SM/WakeFunctionMap.hpp"

vfps::WakeFunctionMap::WakeFunctionMap( std::shared_ptr<PhaseSpace> in
                                      , std::shared_ptr<PhaseSpace> out
                                      , const meshindex_t xsize
                                      , const meshindex_t ysize
                                      , const InterpolationType it
                                      , bool interpol_clamp
                                      , oclhptr_t oclh
                                      )
  : WakeKickMap( in,out,xsize,ysize,it,interpol_clamp,oclh
               #if defined INOVESA_USE_OPENCL and defined INOVESA_USE_OPENGL
               , 0
               #endif // INOVESA_USE_OPENCL and INOVESA_USE_OPENGL
               )
  , _xaxis(Ruler<meshaxis_t>( 2*xsize
                            , in->getMin(0)-in->length(0)/2
                            , in->getMax(0)+in->length(0)/2
                            , in->getScale(0)))
  , _wakefunction(new meshaxis_t[2*xsize])
  , _wakesize(2*xsize)
{
}

vfps::WakeFunctionMap::WakeFunctionMap( std::shared_ptr<PhaseSpace> in
                                      , std::shared_ptr<PhaseSpace> out
                                      , const meshindex_t xsize
                                      , const meshindex_t ysize
                                      , const std::string fname
                                      , const double sigmaE
                                      , const double E0
                                      , const double Ib
                                      , const double dt
                                      , const InterpolationType it
                                      , const bool interpol_clamp
                                      , oclhptr_t oclh
                                      )
  : WakeFunctionMap( in,out,xsize,ysize,it,interpol_clamp,oclh)
{
    /*  1e12*Ib*dt: wake file has charge in pC and for one revolution
     *  1/(ps->getDelta(1)*sigmaE*E0): eV -> pixels
     */
    _wakeFromFile(fname,1e12*Ib*dt/(in->getDelta(1)*E0*sigmaE));
}

vfps::WakeFunctionMap::WakeFunctionMap( std::shared_ptr<PhaseSpace> in
                                      , std::shared_ptr<PhaseSpace> out
                                      , const meshindex_t xsize
                                      , const meshindex_t ysize
                                      , const ElectricField* csr
                                      , const InterpolationType it
                                      , const bool interpol_clamp
                                      , oclhptr_t oclh
                                      )
  : WakeFunctionMap( in,out,xsize,ysize,it,interpol_clamp, oclh)
{
    std::copy_n(csr->getWakefunction(),2*xsize,_wakefunction);
}

vfps::WakeFunctionMap::~WakeFunctionMap() noexcept
{
    delete [] _wakefunction;
    #ifdef INOVESA_ENABLE_CLPROFILING
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
    const integral_t* density = _in->getProjection(0);
    for (unsigned int i=0;i<_xsize;i++) {
        _offset[i] = 0;
        for (unsigned int j=0;j<_xsize;j++) {
            _offset[i] += meshaxis_t(density[j]/charge*_wakefunction[_xsize+i-j]);
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
