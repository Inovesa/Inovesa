#include "HM/WakeFunctionMap.hpp"

vfps::WakeFunctionMap::WakeFunctionMap(vfps::PhaseSpace* in, vfps::PhaseSpace* out,
                const vfps::meshindex_t xsize, const vfps::meshindex_t ysize,
                const vfps::HeritageMap::InterpolationType it) :
    WakeKickMap(in,out,xsize,ysize,it),
    _wakefunction(new meshaxis_t[2*xsize]),
    _wakesize(2*xsize)
{
}

vfps::WakeFunctionMap::WakeFunctionMap(vfps::PhaseSpace* in, vfps::PhaseSpace* out,
                const vfps::meshindex_t xsize, const vfps::meshindex_t ysize,
                const std::vector<std::pair<meshaxis_t, double>> wakefunction,
                const InterpolationType it) :
    WakeFunctionMap(in,out,xsize,ysize,it)
{
    const Ruler<meshaxis_t> xaxis(_wakesize,2*in->getMin(0),2*in->getMax(0));
    size_t x_other=1;
    size_t x_mine=0;
    std::array<interpol_t,4> qic;
    bool smallwake=false;

    if (xaxis[0] < wakefunction[1].first) {
        smallwake = true;
    }

    while (x_mine < _wakesize) {
        while (xaxis[x_mine] > wakefunction[x_other+1].first) {
            if (x_other+3 < wakefunction.size()) {
                x_other++;
            } else {
                smallwake = true;
                break;
            }
        }
        calcCoefficiants(qic.data(),
                         xaxis[x_mine]-wakefunction[x_other].first,4);
        _wakefunction[x_mine]   = qic[0]*wakefunction[x_other-1].second
                                + qic[1]*wakefunction[x_other  ].second
                                + qic[2]*wakefunction[x_other+1].second
                                + qic[3]*wakefunction[x_other+2].second;
        x_mine++;
    }
    if (smallwake) {
        Display::printText("Warning: Given wake to small.");
    }
}

vfps::WakeFunctionMap::WakeFunctionMap(vfps::PhaseSpace* in, vfps::PhaseSpace* out,
                const vfps::meshindex_t xsize, const vfps::meshindex_t ysize,
                               const vfps::ElectricField* csr,
                               const vfps::HeritageMap::InterpolationType it) :
    WakeFunctionMap(in,out,xsize,ysize,it)
{
    std::copy_n(csr->getWakefunction(),2*xsize,_wakefunction);
}

vfps::WakeFunctionMap::~WakeFunctionMap()
{
    delete [] _wakefunction;
}

void vfps::WakeFunctionMap::update()
{
    #if INOVESA_USE_CL
    if (OCLH::active) {
    _in->syncCLMem(clCopyDirection::dev2cpu);
    }
    #endif
    integral_t charge = _in->integral();
    const integral_t* density = _in->getProjection(0);
    for (unsigned int i=0;i<_xsize;i++) {
        _force[i] = 0;
        for (unsigned int j=0;j<_xsize;j++) {
            _force[i] += meshaxis_t(density[j]/charge*_wakefunction[_xsize+i-j]);
        }
    }
}
