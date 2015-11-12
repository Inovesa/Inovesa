#include "HM/WakePotentialMap.hpp"

vfps::WakePotentialMap::WakePotentialMap(
        vfps::PhaseSpace *in, vfps::PhaseSpace *out,
        const vfps::meshindex_t xsize, const vfps::meshindex_t ysize,
        ElectricField *field, const vfps::HeritageMap::InterpolationType it):
    WakeKickMap(in,out,xsize,ysize,it),
    _field(field)
{
}

void vfps::WakePotentialMap::update()
{
    #if INOVESA_USE_CL
    if (OCLH::active) {
    _in->syncCLMem(PhaseSpace::clCopyDirection::dev2cpu);
    }
    #endif
    std::copy_n(_field->wakePotential(),_xsize,_force);
    for (meshindex_t x=0;x<_xsize;x++) {
        _force[x] *= 2e-8;
    }
}
