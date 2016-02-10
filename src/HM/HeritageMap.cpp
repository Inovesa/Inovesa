/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   *
 * Copyright (c) 2014-2015: Patrik Schönfeldt                                 *
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

#include "HM/HeritageMap.hpp"

vfps::HeritageMap::HeritageMap(PhaseSpace* in, PhaseSpace* out,
                               meshindex_t xsize, meshindex_t ysize,
                               size_t memsize,
                               uint_fast8_t interpoints,
                               uint_fast8_t intertype) :
    _ip(interpoints),
    _it(intertype),
    _hinfo(new hi[std::max(memsize,static_cast<size_t>(16))]),
    _size(xsize*ysize),
    _xsize(xsize),
    _ysize(ysize),
    _in(in),
    _out(out)
{
    #ifdef INOVESA_USE_CL
    _cl_code  += "typedef struct { uint src; data_t weight; } hi;\n";
    #endif // INOVESA_USE_CL
}

vfps::HeritageMap::HeritageMap(PhaseSpace* in, PhaseSpace* out,
                               size_t xsize, size_t ysize,
                               uint_fast8_t interpoints,
                               uint_fast8_t intertype) :
    HeritageMap(in,out,xsize,ysize,xsize*ysize*interpoints,
                interpoints,intertype)
{
}

vfps::HeritageMap::~HeritageMap()
{
    delete [] _hinfo;
}

void vfps::HeritageMap::apply()
{
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        #ifdef INOVESA_SYNC_CL
        _in->syncCLMem(clCopyDirection::cpu2dev);
        #endif // INOVESA_SYNC_CL
        OCLH::queue.enqueueNDRangeKernel (
                    applyHM,
                    cl::NullRange,
                    cl::NDRange(_size));
        #ifdef CL_VERSION_1_2
        OCLH::queue.enqueueBarrierWithWaitList();
        #else // CL_VERSION_1_2
        OCLH::queue.enqueueBarrier();
        #endif // CL_VERSION_1_2
        #ifdef INOVESA_SYNC_CL
        _out->syncCLMem(clCopyDirection::dev2cpu);
        #endif // INOVESA_SYNC_CL
    } else
    #endif // INOVESA_USE_CL
    {
        meshdata_t* data_in = _in->getData();
        meshdata_t* data_out = _out->getData();

        for (meshindex_t i=0; i< _size; i++) {
            data_out[i] = 0;
            for (meshindex_t j=0; j<_ip; j++) {
                hi h = _hinfo[i*_ip+j];
                data_out[i] += data_in[h.index]*static_cast<meshdata_t>(h.weight);
            }
        }
    }
}


#ifdef INOVESA_USE_CL
void vfps::HeritageMap::genCode4HM1D()
{
    _cl_code += R"(
    __kernel void applyHM1D(const __global data_t* src,
                            const __global hi* hm,
                            const uint hm_len,
                            __global data_t* dst)
    {
        data_t value = 0;
        const uint i = get_global_id(0);
        const uint offset = i*hm_len;
        for (uint j=0; j<hm_len; j++)
        {
            value += mult(src[hm[offset+j].src],hm[offset+j].weight);
        }
        dst[i] = value;
    }
)";
}
#endif // INOVESA_USE_CL

void vfps::HeritageMap::calcCoefficiants(vfps::interpol_t* ic,
                                         const vfps::interpol_t f,
                                         const uint_fast8_t it) const
{
    switch (it) {
    case InterpolationType::none:
        ic[0] = 1;
        break;
    case InterpolationType::linear:
        ic[0] = interpol_t(1)-f;
        ic[1] = f;
        break;
    case InterpolationType::quadratic:
        ic[0] = f*(f-interpol_t(1))/interpol_t(2);
        ic[1] = interpol_t(1)-f*f;
        ic[2] = f*(f+interpol_t(1))/interpol_t(2);
        break;
    case InterpolationType::cubic:
        ic[0] = (f-interpol_t(1))*(f-interpol_t(2))*f
                * interpol_t(-1./6.);
        ic[1] = (f+interpol_t(1))*(f-interpol_t(1))
                * (f-interpol_t(2)) / interpol_t(2);
        ic[2] = (interpol_t(2)-f)*f*(f+interpol_t(1))
                / interpol_t(2);
        ic[3] = f*(f+interpol_t(1))*(f-interpol_t(1))
                * interpol_t(1./6.);
        break;
    }
}

void vfps::HeritageMap::notBoundMessage()
{
    vfps::Display::printText("\tBound interpolation not implemented for this scheme.");
}
