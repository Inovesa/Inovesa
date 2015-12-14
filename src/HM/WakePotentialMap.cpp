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
    _field->wakePotential();
    OCLH::queue.enqueueCopyBuffer(_field->_wakepotential_buf,_force_buf,
                                  0,0,sizeof(meshaxis_t)*_xsize);
    } else
    #endif
    {
    std::copy_n(_field->wakePotential(),_xsize,_force.data());
    }
}
