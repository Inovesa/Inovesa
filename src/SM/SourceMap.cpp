// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * This file is part of Inovesa (github.com/Inovesa/Inovesa).
 * It's copyrighted by the contributors recorded
 * in the version control history of the file.
 */

#include "SM/SourceMap.hpp"

#include "MessageStrings.hpp"

vfps::SourceMap::SourceMap( std::shared_ptr<PhaseSpace> in
                          , std::shared_ptr<PhaseSpace> out
                          , meshindex_t xsize
                          , meshindex_t ysize
                          , size_t memsize
                          , uint_fast8_t interpoints
                          , uint_fast8_t intertype
                          , oclhptr_t oclh
                          )
  : _ip(interpoints)
  , _it(intertype)
  , _hinfo(new hi[std::max(memsize,static_cast<size_t>(16))])
  , _size(xsize*ysize)
  , _xsize(xsize)
  , _ysize(ysize)
  #if INOVESA_ENABLE_CLPROFILING == 1
  , applySMEvents(std::make_unique<std::vector<cl::Event*>>())
  , syncSMEvents(std::make_unique<std::vector<cl::Event*>>())
  #endif // INOVESA_ENABLE_CLPROFILING
  , _axis(std::array<meshRuler_ptr,2>{{in->getAxis(0),in->getAxis(1)}})
  , _in(in)
  , _out(out)
  , _oclh(oclh)
{
    #if INOVESA_USE_OPENCL == 1
    _cl_code  += "typedef struct { uint src; data_t weight; } hi;\n";
    #endif // INOVESA_USE_OPENCL
}

vfps::SourceMap::SourceMap( std::shared_ptr<PhaseSpace> in
                          , std::shared_ptr<PhaseSpace> out
                          , size_t xsize, size_t ysize
                          , uint_fast8_t interpoints
                          , uint_fast8_t intertype
                          , oclhptr_t oclh
                          )
  : SourceMap( in,out,xsize,ysize,xsize*ysize*interpoints
             , interpoints,intertype,oclh
             )
{
}

vfps::SourceMap::~SourceMap() noexcept
{
    delete [] _hinfo;
    #if INOVESA_ENABLE_CLPROFILING == 1
    for (auto ev : *applySMEvents) {
        delete ev;
    }
    for (auto ev : *syncSMEvents) {
        delete ev;
    }
    #endif // INOVESA_ENABLE_CLPROFILING
}


#if INOVESA_ENABLE_CLPROFILING == 1
void vfps::SourceMap::saveTimings(std::string mapname) {
    if (_oclh) {
        _oclh->saveTimings(applySMEvents.get(),"Apply"+mapname);
        _oclh->saveTimings(syncSMEvents.get(),"Sync"+mapname);
    }
}
#endif // INOVESA_ENABLE_CLPROFILING

void vfps::SourceMap::apply()
{
    #if INOVESA_USE_OPENCL == 1
    if (_oclh) {
        #if INOVESA_SYNC_CL == 1
        _in->syncCLMem(OCLH::clCopyDirection::cpu2dev);
        #endif // INOVESA_SYNC_CL
        _oclh->enqueueNDRangeKernel( applySM
                                  , cl::NullRange
                                  , cl::NDRange(_size)
                                  #if INOVESA_ENABLE_CLPROFILING == 1
                                  , cl::NullRange
                                  , nullptr
                                  , nullptr
                                  , applySMEvents.get()
                                  #endif // INOVESA_ENABLE_CLPROFILING
                                  );
        #if INOVESA_SYNC_CL == 1
        _out->syncCLMem(OCLH::clCopyDirection::dev2cpu);
        #endif // INOVESA_SYNC_CL
    } else
    #endif // INOVESA_USE_OPENCL
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

void vfps::SourceMap::applyTo(std::vector<vfps::PhaseSpace::Position> &particles)
{
    for (PhaseSpace::Position& particle : particles) {
        particle = apply(particle);
    }
}


#if INOVESA_USE_OPENCL == 1
void vfps::SourceMap::genCode4SM1D()
{
    _cl_code += R"(
    __kernel void applySM1D(const __global data_t* src,
                            const __global hi* sm,
                            const uint sm_len,
                            __global data_t* dst)
    {
        data_t value = 0;
        const uint i = get_global_id(0);
        const uint offset = i*sm_len;
        for (uint j=0; j<sm_len; j++)
        {
            value += mult(src[sm[offset+j].src],sm[offset+j].weight);
        }
        dst[i] = value;
    }
)";
}
#endif // INOVESA_USE_OPENCL

void vfps::SourceMap::calcCoefficiants(vfps::interpol_t* ic,
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

void vfps::SourceMap::notClampedMessage()
{
    Display::printText("WARNING: Clamped interpolation not implemented "
                       "for this interpolation scheme.");
}
