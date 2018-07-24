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

#include "SM/FokkerPlanckMap.hpp"

vfps::FokkerPlanckMap::FokkerPlanckMap( std::shared_ptr<PhaseSpace> in
                                      , std::shared_ptr<PhaseSpace> out
                                      , const meshindex_t xsize
                                      , const meshindex_t ysize
                                      , FPType fptype, FPTracking fptrack
                                      , timeaxis_t e1
                                      , DerivationType dt
                                      , oclhptr_t oclh
                                      )
  : SourceMap( in, out, 1, ysize, dt, dt, oclh)
  , _dampdecr(e1)
  , _prng(std::mt19937(std::random_device{}()))
  , _normdist( std::normal_distribution<meshaxis_t>( 0
                                                   , std::sqrt(2*e1)
                                                     / in->getDelta(1)))
  , _fptrack(fptrack)
  , _fptype(fptype)
  , _meshxsize(xsize)
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

            if (_fptype != FPType::none && _fptype != FPType::diffusion_only) {
                const meshaxis_t pos = in->p(j);
                _hinfo[j*_ip  ].weight += -e1_2d*pos;
                _hinfo[j*_ip+1].weight +=  e1;
                _hinfo[j*_ip+2].weight +=  e1_2d*pos;
            }
            if (_fptype != FPType::none && _fptype != FPType::damping_only) {
                _hinfo[j*_ip  ].weight +=                e1_d2;
                _hinfo[j*_ip+1].weight += interpol_t(-2)*e1_d2;
                _hinfo[j*_ip+2].weight +=                e1_d2;
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
            const meshaxis_t pos = in->p(j);
            _hinfo[j*_ip  ]={j-2,0};
            _hinfo[j*_ip+1]={j-1,0};
            _hinfo[j*_ip+2]={j  ,1};
            _hinfo[j*_ip+3]={j+1,0};
            if (_fptype != FPType::none && _fptype != FPType::diffusion_only) {
                _hinfo[j*_ip  ].weight +=    e1_6d*interpol_t( 1)*pos;
                _hinfo[j*_ip+1].weight +=    e1_6d*interpol_t(-6)*pos;
                _hinfo[j*_ip+2].weight += e1+e1_6d*interpol_t( 3)*pos;
                _hinfo[j*_ip+3].weight +=    e1_6d*interpol_t( 2)*pos;
            }
            if (_fptype != FPType::none && _fptype != FPType::damping_only) {
                _hinfo[j*_ip+1].weight +=    e1_d2;
                _hinfo[j*_ip+2].weight += interpol_t(-2)*e1_d2;
                _hinfo[j*_ip+3].weight +=    e1_d2;
            }
        }
        for (meshindex_t j=ycenter; j<static_cast<meshindex_t>(_ysize-2);j++) {
            const meshaxis_t pos = in->p(j);
            _hinfo[j*_ip  ]={j-1,0};
            _hinfo[j*_ip+1]={j  ,1};
            _hinfo[j*_ip+2]={j+1,0};
            _hinfo[j*_ip+3]={j+2,0};

            if (_fptype != FPType::none && _fptype != FPType::diffusion_only) {
                _hinfo[j*_ip  ].weight +=    e1_6d*interpol_t(-2)*pos;
                _hinfo[j*_ip+1].weight += e1+e1_6d*interpol_t(-3)*pos;
                _hinfo[j*_ip+2].weight +=    e1_6d*interpol_t( 6)*pos;
                _hinfo[j*_ip+3].weight +=    e1_6d*interpol_t(-1)*pos;
            }
            if (_fptype != FPType::none && _fptype != FPType::damping_only) {
                _hinfo[j*_ip  ].weight +=                e1_d2;
                _hinfo[j*_ip+1].weight += interpol_t(-2)*e1_d2;
                _hinfo[j*_ip+2].weight +=                e1_d2;
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

    #ifdef INOVESA_USE_OPENCL
    if (_oclh) {
    _cl_code += R"(
    __kernel void applySM_Y(const __global data_t* src,
                            const __global hi* sm,
                            const uint sm_len,
                            const uint ysize,
                            __global data_t* dst)
    {
        data_t value = 0;
        const uint x = get_global_id(0);
        const uint y = get_global_id(1);
        const uint smoffset = y*sm_len;
        const uint meshoffs = x*ysize;
        for (uint j=0; j<sm_len; j++)
        {
            value += mult(  src[meshoffs+sm[smoffset+j].src],
                            sm[smoffset+j].weight);
        }
        dst[meshoffs+y] = value;
    }
    )";

    _cl_prog = _oclh->prepareCLProg(_cl_code);

    if (_oclh) {
        _sm_buf = cl::Buffer(_oclh->context,
                             CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                             sizeof(hi)*_ip*_ysize,
                             _hinfo);
        applySM = cl::Kernel(_cl_prog, "applySM_Y");
        applySM.setArg(0, _in->data_buf);
        applySM.setArg(1, _sm_buf);
        applySM.setArg(2, _ip);
        applySM.setArg(3, _ysize);
        applySM.setArg(4, _out->data_buf);
    }
    }
#endif
}

vfps::FokkerPlanckMap::~FokkerPlanckMap() noexcept
#ifdef INOVESA_ENABLE_CLPROFILING
{
    saveTimings("FokkerPlanckMap");
}
#else
    = default;
#endif // INOVESA_ENABLE_CLPROFILING

void vfps::FokkerPlanckMap::apply()
{
    #ifdef INOVESA_USE_OPENCL
    if (_oclh) {
        #ifdef INOVESA_SYNC_CL
        _in->syncCLMem(OCLH::clCopyDirection::cpu2dev);
        #endif // INOVESA_SYNC_CL
        _oclh->enqueueNDRangeKernel( applySM
                                  , cl::NullRange
                                  , cl::NDRange(_meshxsize,_ysize)
                                  #ifdef INOVESA_ENABLE_CLPROFILING
                                  , cl::NullRange
                                  , nullptr
                                  , nullptr
                                  , applySMEvents.get()
                                  #endif // INOVESA_ENABLE_CLPROFILING
                                  );
        _oclh->enqueueBarrier();
        #ifdef INOVESA_SYNC_CL
        _out->syncCLMem(OCLH::clCopyDirection::dev2cpu);
        #endif // INOVESA_SYNC_CL
    } else
    #endif // INOVESA_USE_OPENCL
    {
        meshdata_t* data_in = _in->getData();
        meshdata_t* data_out = _out->getData();

        for (meshindex_t x=0; x< _meshxsize; x++) {
            const meshindex_t offs = x*_ysize;
            for (meshindex_t y=0; y< _ysize; y++) {
                meshdata_t value = 0;
                for (uint_fast8_t j=0; j<_ip; j++) {
                    hi h = _hinfo[y*_ip+j];
                    value += data_in[offs+h.index]
                                     *  static_cast<meshdata_t>(h.weight);
                }
                data_out[offs+y] = value;
            }
        }
    }
}

vfps::PhaseSpace::Position
vfps::FokkerPlanckMap::apply(PhaseSpace::Position pos) const
{
    switch (_fptrack){
    case FPTracking::none:
    default:
        break;
    case FPTracking::approximation1:
        {
        meshindex_t yi = std::min( static_cast<meshindex_t>(std::floor(pos.y))
                                 , _ysize);
        interpol_t offset = 0;

        for (uint_fast8_t j=0; j<_ip; j++) {
            hi h = _hinfo[yi*_ip+j];
            interpol_t dy = static_cast<interpol_t>(yi)
                          - static_cast<interpol_t>(h.index);
            offset += dy*h.weight;
        }
        pos.y = std::max(static_cast<meshaxis_t>(1),
        std::min(pos.y+offset,static_cast<meshaxis_t>(_ysize-1)));
        }
        break;
    case FPTracking::approximation2:
        {
        meshdata_t* data_in = _in->getData();
        std::make_signed<meshindex_t>::type xi
            = std::min( static_cast<decltype(_in->nMeshCells(0))>(std::floor(pos.x))
                      , _in->nMeshCells(0)-1);
        std::make_signed<meshindex_t>::type yi
            = std::min( static_cast<decltype(_ysize)>(std::floor(pos.y))
                      , _ysize-1);
        const meshindex_t offs = xi*_ysize;
        interpol_t offset = 0;
        meshdata_t charge = 0;

        for (uint_fast8_t j=0; j<_ip; j++) {
            hi h = _hinfo[yi*_ip+j];
            charge += data_in[offs+h.index]*h.weight;
            offset += data_in[offs+h.index]*h.weight
                    * (static_cast<std::make_signed<meshindex_t>::type>(h.index)
                      - yi);
        }
        offset /= charge;
        pos.y = std::max( static_cast<meshaxis_t>(1)
                        , std::min(pos.y+offset
                                  , static_cast<meshaxis_t>(_ysize-1)));
        }
        break;
    case FPTracking::stochastic:
        auto delta = pos.y*_dampdecr+_normdist(_prng);
        pos.y -= delta;
        break;
    }
    return pos;
}

