#include "HM/WakeImpedanceMap.hpp"




vfps::WakeImpedanceMap::WakeImpedanceMap(vfps::PhaseSpace *in, vfps::PhaseSpace *out,
        const vfps::meshindex_t xsize, const vfps::meshindex_t ysize,
        ElectricField *field, const vfps::Impedance* impedance,
        const vfps::HeritageMap::InterpolationType it):
    WakeKickMap(in,out,xsize,ysize,it),
    _field(field),
    _impedance(impedance)
{
}

void vfps::WakeImpedanceMap::apply()
{
    #if INOVESA_USE_CL
    if (OCLH::active) {
    _in->syncCLMem(PhaseSpace::clCopyDirection::dev2cpu);
    }
    #endif
    integral_t charge = _in->integral();
    // get wake potential from ElectricField
    meshaxis_t* wake = _field->wakePotential();

    for (unsigned int i=0;i<_xsize;i++) {
        _force[i] += wake[i]/charge; // renormalize to units of p (2DO)
    }
    WakeKickMap::apply();
}
