#include "HM/WakePotentialMap.hpp"

vfps::WakePotentialMap::WakePotentialMap(
        vfps::PhaseSpace *in, vfps::PhaseSpace *out,
        const vfps::meshindex_t xsize, const vfps::meshindex_t ysize,
        ElectricField *field, const vfps::HeritageMap::InterpolationType it):
    WakeKickMap(in,out,xsize,ysize,it),
    _field(field)
{
    #ifdef INOVESA_USE_CL
    _force_buf = _field->_wakepotential_buf;
    #endif
}

void vfps::WakePotentialMap::update()
{
    #if INOVESA_USE_CL
    if (OCLH::active) {
    _field->wakePotential();
    } else
    #endif
    {
    std::copy_n(_field->wakePotential(),_xsize,_force.data());
    }
}
