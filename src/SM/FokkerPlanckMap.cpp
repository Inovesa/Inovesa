/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
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

#include "SM/FokkerPlanckMap.hpp"

vfps::FokkerPlanckMap::FokkerPlanckMap(PhaseSpace* in, PhaseSpace* out,
                                       const meshindex_t xsize,
                                       const meshindex_t ysize,
                                       FPType fpt, timeaxis_t e1,
                                       DerivationType dt)
    :
    SourceMap(in, out, 1, ysize, dt, dt),
    _meshxsize(xsize)
{
    // the following doubles should be interpol_t
    const interpol_t e1_2d = e1/(interpol_t(2)*in->getDelta(1));
    const interpol_t e1_6d = e1/(interpol_t(6)*in->getDelta(1));
    const interpol_t e1_d2 = e1/(in->getDelta(1)*in->getDelta(1));

    const meshaxis_t ycenter = in->getAxis(1)->zerobin();

    switch (dt) {
    case DerivationType::two_sided:
        _hinfo[0] = {0,0};
        _hinfo[1] = {0,0};
        _hinfo[2] = {0,0};
        for (meshindex_t j=1; j< _ysize-1; j++) {
            _hinfo[j*_ip  ]={j-1,0};
            _hinfo[j*_ip+1]={j  ,1};
            _hinfo[j*_ip+2]={j+1,0};

            if (fpt == FPType::full || fpt == FPType::damping_only) {
                const meshaxis_t pos = in->x(1,j);
                _hinfo[j*_ip  ].weight += -e1_2d*pos;
                _hinfo[j*_ip+1].weight +=  e1;
                _hinfo[j*_ip+2].weight +=  e1_2d*pos;
            }
            if (fpt == FPType::full || fpt == FPType::diffusion_only) {
                _hinfo[j*_ip  ].weight +=                 e1_d2;
                _hinfo[j*_ip+1].weight += interpol_t(-2)*e1_d2;
                _hinfo[j*_ip+2].weight +=                 e1_d2;
            }
        }
        _hinfo[(_ysize-1)*_ip+0] = {0,0};
        _hinfo[(_ysize-1)*_ip+1] = {0,0};
        _hinfo[(_ysize-1)*_ip+2] = {0,0};
        break;
    case DerivationType::cubic:
        _hinfo[0] = {0,0};
        _hinfo[1] = {0,0};
        _hinfo[2] = {0,0};
        _hinfo[3] = {0,0};
        _hinfo[_ip+0] = {0,0};
        _hinfo[_ip+1] = {0,0};
        _hinfo[_ip+2] = {0,0};
        _hinfo[_ip+3] = {0,0};
        for (meshindex_t j=2; j< ycenter; j++) {
            const meshaxis_t pos = in->x(1,j);
            _hinfo[j*_ip  ]={j-2,0};
            _hinfo[j*_ip+1]={j-1,0};
            _hinfo[j*_ip+2]={j  ,1};
            _hinfo[j*_ip+3]={j+1,0};
            if (fpt == FPType::full || fpt == FPType::damping_only) {
                _hinfo[j*_ip  ].weight +=    e1_6d*interpol_t( 1)*pos;
                _hinfo[j*_ip+1].weight +=    e1_6d*interpol_t(-6)*pos;
                _hinfo[j*_ip+2].weight += e1+e1_6d*interpol_t( 3)*pos;
                _hinfo[j*_ip+3].weight +=    e1_6d*interpol_t( 2)*pos;
            }
            if (fpt == FPType::full || fpt == FPType::diffusion_only) {
                _hinfo[j*_ip+1].weight +=    e1_d2;
                _hinfo[j*_ip+2].weight += interpol_t(-2)*e1_d2;
                _hinfo[j*_ip+3].weight +=    e1_d2;
            }
        }
        for (meshindex_t j=ycenter; j<static_cast<meshindex_t>(_ysize-2);j++) {
            const meshaxis_t pos = in->x(1,j);
            _hinfo[j*_ip  ]={j-1,0};
            _hinfo[j*_ip+1]={j  ,1};
            _hinfo[j*_ip+2]={j+1,0};
            _hinfo[j*_ip+3]={j+2,0};

            if (fpt == FPType::full || fpt == FPType::damping_only) {
                _hinfo[j*_ip  ].weight +=    e1_6d*interpol_t(-2)*pos;
                _hinfo[j*_ip+1].weight += e1+e1_6d*interpol_t(-3)*pos;
                _hinfo[j*_ip+2].weight +=    e1_6d*interpol_t( 6)*pos;
                _hinfo[j*_ip+3].weight +=    e1_6d*interpol_t(-1)*pos;
            }
            if (fpt == FPType::full || fpt == FPType::diffusion_only) {
                _hinfo[j*_ip  ].weight +=                 e1_d2;
                _hinfo[j*_ip+1].weight += interpol_t(-2)*e1_d2;
                _hinfo[j*_ip+2].weight +=                 e1_d2;
            }
        }
        _hinfo[(_ysize-2)*_ip  ] = {0,0};
        _hinfo[(_ysize-2)*_ip+1] = {0,0};
        _hinfo[(_ysize-2)*_ip+2] = {0,0};
        _hinfo[(_ysize-2)*_ip+3] = {0,0};
        _hinfo[(_ysize-1)*_ip+0] = {0,0};
        _hinfo[(_ysize-1)*_ip+1] = {0,0};
        _hinfo[(_ysize-1)*_ip+2] = {0,0};
        _hinfo[(_ysize-1)*_ip+3] = {0,0};
        break;
    }

    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
    _cl_code += R"(
    __kernel void applyHM_Y(const __global data_t* src,
                            const __global hi* hm,
                            const uint hm_len,
                            const uint ysize,
                            __global data_t* dst)
    {
        data_t value = 0;
        const uint x = get_global_id(0);
        const uint y = get_global_id(1);
        const uint hmoffset = y*hm_len;
        const uint meshoffs = x*ysize;
        for (uint j=0; j<hm_len; j++)
        {
            value += mult(  src[meshoffs+hm[hmoffset+j].src],
                            hm[hmoffset+j].weight);
        }
        dst[meshoffs+y] = value;
    }
    )";

    _cl_prog = OCLH::prepareCLProg(_cl_code);

    if (OCLH::active) {
        _hi_buf = cl::Buffer(OCLH::context,
                             CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                             sizeof(hi)*_ip*_ysize,
                             _hinfo);
        applyHM = cl::Kernel(_cl_prog, "applyHM_Y");
        applyHM.setArg(0, _in->data_buf);
        applyHM.setArg(1, _hi_buf);
        applyHM.setArg(2, _ip);
        applyHM.setArg(3, _ysize);
        applyHM.setArg(4, _out->data_buf);
    }
    }
    #endif
}

void vfps::FokkerPlanckMap::apply()
{
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        #ifdef INOVESA_SYNC_CL
        _in->syncCLMem(clCopyDirection::cpu2dev);
        #endif // INOVESA_SYNC_CL
        OCLH::queue.enqueueNDRangeKernel (
                    applyHM,
                    cl::NullRange,
                    cl::NDRange(_meshxsize,_ysize));
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

        for (meshindex_t x=0; x< _meshxsize; x++) {
            const meshindex_t offs = x*_ysize;
            for (meshindex_t y=0; y< _ysize; y++) {
                data_out[offs+y] = 0;
                for (uint_fast8_t j=0; j<_ip; j++) {
                    hi h = _hinfo[y*_ip+j];
                    data_out[offs+y] += data_in[offs+h.index]
                                     *  static_cast<meshdata_t>(h.weight);
                }
            }
        }
    }
}

vfps::PhaseSpace::position
vfps::FokkerPlanckMap::apply(PhaseSpace::position pos) const
{
    meshindex_t yi = std::min(static_cast<meshindex_t>(std::floor(pos.y)),_ysize);
    interpol_t offset = 0;

    for (uint_fast8_t j=0; j<_ip; j++) {
        hi h = _hinfo[yi*_ip+j];
        interpol_t dy = static_cast<interpol_t>(yi)
                      - static_cast<interpol_t>(h.index);
        offset += dy*h.weight;
    }
    pos.y = std::max(static_cast<meshaxis_t>(1),
                     std::min(pos.y+offset,static_cast<meshaxis_t>(_ysize-1)));
    return pos;
}

