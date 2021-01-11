// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#include <stdexcept>
#include <type_traits>

#include "SM/RotationMap.hpp"

vfps::RotationMap::RotationMap(std::shared_ptr<PhaseSpace> in,
                               std::shared_ptr<PhaseSpace> out,
                               const meshindex_t xsize,
                               const meshindex_t ysize,
                               const meshaxis_t angle,
                               const InterpolationType it,
                               const bool interpol_clamped,
                               const meshindex_t rotmapsize,
                               oclhptr_t oclh) :
    SourceMap( in,out,xsize,ysize,size_t(rotmapsize)*it*it,it*it,it,oclh),
    _rotmapsize(rotmapsize),
    _clamp(interpol_clamped),
    _cos_dt(std::cos(-angle)),
    _sin_dt(std::sin(-angle))
{
    if (_rotmapsize == 0) {
        #if INOVESA_USE_OPENCL == 1
        if (_oclh) {
            rot = {{float(_cos_dt),float(_sin_dt)}};
            imgsize = {{cl_int(_xsize),cl_int(_ysize)}};
            zerobin = {{(_axis[0]->zerobin()),(_axis[1]->zerobin())}};


            genCode4Rotation();
            _cl_prog  = _oclh->prepareCLProg(_cl_code);

            applySM = cl::Kernel(_cl_prog, "applyRotation");
            applySM.setArg(0, _in->data_buf);
            applySM.setArg(1, rot);
            applySM.setArg(2, imgsize);
            applySM.setArg(3, zerobin);
            applySM.setArg(4, _out->data_buf);
        }
        #endif // INOVESA_USE_OPENCL
    } else {
        for (meshindex_t q_i=0; q_i< _xsize; q_i++) {
            for(meshindex_t p_i=0; p_i< _ysize; p_i++) {
                genHInfo(q_i,p_i,&_hinfo[(q_i*ysize+p_i)*_ip]);
            }
        }
        #if INOVESA_USE_OPENCL == 1
        if (_oclh) {
            _sm_buf = cl::Buffer(_oclh->context,
                                 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                 sizeof(hi)*_ip*_rotmapsize,
                                 _hinfo.data());
            if (_clamp) {
                if (it == InterpolationType::cubic) {
                    if (_rotmapsize == PhaseSpace::nxy) {
                        genCode4SM4sat();
                        _cl_prog  = _oclh->prepareCLProg(_cl_code);

                        applySM = cl::Kernel(_cl_prog, "applySM4sat");
                        applySM.setArg(0, _in->data_buf);
                        applySM.setArg(1, _sm_buf);
                        applySM.setArg(2, _out->data_buf);
                    }
                }
            } else {
                genCode4SM1D();
                _cl_prog  = _oclh->prepareCLProg(_cl_code);

                applySM = cl::Kernel(_cl_prog, "applySM1D");
                applySM.setArg(0, _in->data_buf);
                applySM.setArg(1, _sm_buf);
                applySM.setArg(2, _ip);
                applySM.setArg(3, _out->data_buf);
            }
        }
        #endif // INOVESA_USE_OPENCL
    }
    if (interpol_clamped && !(it==InterpolationType::cubic && rotmapsize>0)) {
        throw std::invalid_argument("Clamping only supported"
                                    "with cubic interpolation"
                                    "and rotmapsize > 0.");
    }
}

#if INOVESA_ENABLE_CLPROFILING == 1
vfps::RotationMap::~RotationMap() noexcept
{
    saveTimings("RotationMap");
}
#endif // INOVESA_ENABLE_CLPROFILING

void vfps::RotationMap::apply()
{
    #if INOVESA_USE_OPENCL == 1
    if (_oclh) {
        #if INOVESA_SYNC_CL == 1
        _in->syncCLMem(clCopyDirection::cpu2dev);
        #endif // INOVESA_SYNC_CL
        if (_rotmapsize == 0) {
             // stay away from mesh borders
            _oclh->enqueueNDRangeKernel( applySM
                                      , cl::NDRange(1,1)
                                      , cl::NDRange(_xsize-_it+1,_ysize-_it+1)
                                      #if INOVESA_ENABLE_CLPROFILING == 1
                                      , cl::NullRange
                                      , nullptr
                                      , nullptr
                                      , applySMEvents.get()
                                      #endif // INOVESA_ENABLE_CLPROFILING
                                      );
        } else {
            _oclh->enqueueNDRangeKernel( applySM
                                      , cl::NullRange
                                      , cl::NDRange(_rotmapsize)
                                      #if INOVESA_ENABLE_CLPROFILING == 1
                                      , cl::NullRange
                                      , nullptr
                                      ,nullptr
                                      , applySMEvents.get()
                                      #endif // INOVESA_ENABLE_CLPROFILING
                                      );
        }
        _oclh->enqueueBarrier();
        #ifdef INOVESA_SYNC_CL
        _out->syncCLMem(clCopyDirection::dev2cpu);
        #endif // INOVESA_SYNC_CL
    } else
    #endif // INOVESA_USE_OPENCL
    {
        meshdata_t* data_in = _in->getData();
        meshdata_t* data_out = _out->getData();

        if (_rotmapsize == 0) {
            for (meshindex_t q_i=0; q_i< _xsize; q_i++) {
                for(meshindex_t p_i=0; p_i< _ysize; p_i++) {
                    meshindex_t i = q_i*_ysize+p_i;
                    data_out[i] = 0;
                    genHInfo(q_i,p_i,&(_hinfo[0]));
                    for (std::remove_const<decltype(_ip)>::type j=0; j<_ip; j++) {
                        hi h = _hinfo[j];
                        data_out[i] += data_in[h.index]
                                    * static_cast<meshdata_t>(h.weight);
                    }
                }
            }
        } else {
            for (meshindex_t i=0; i< _rotmapsize; i++) {
                data_out[i] = 0;
                for (std::remove_const<decltype(_ip)>::type j=0; j<_ip; j++) {
                    hi h = _hinfo[i*_ip+j];
                    data_out[i] += data_in[h.index]
                                * static_cast<meshdata_t>(h.weight);
                }
                if (_clamp) {
                    // handle overshooting
                    meshdata_t ceil=std::numeric_limits<meshdata_t>::min();
                    meshdata_t flor=std::numeric_limits<meshdata_t>::max();
                    for (meshindex_t x=1; x<=2; x++) {
                        for (meshindex_t y=1; y<=2; y++) {
                            ceil = std::max(
                                ceil, data_in[_hinfo[i*_ip+x*_it+y].index]);
                            flor = std::min(
                                flor, data_in[_hinfo[i*_ip+x*_it+y].index]);
                        }
                    }
                    data_out[i] = std::max(std::min(ceil,data_out[i]),flor);
                }
            }
        }
    }
}

void vfps::RotationMap::applyTo(PhaseSpace::Position& pos) const
{
    PhaseSpace::Position tmp;
    tmp.x = _cos_dt*meshaxis_t(pos.x-_axis[0]->zerobin())
         + _sin_dt*meshaxis_t(pos.y-_axis[1]->zerobin())
         + _axis[0]->zerobin();
    tmp.y = _cos_dt*meshaxis_t(pos.y-_axis[1]->zerobin())
         - _sin_dt*meshaxis_t(pos.x-_axis[0]->zerobin())
         + _axis[1]->zerobin();
    pos = tmp;
}

void vfps::RotationMap::genHInfo(vfps::meshindex_t x0,
                                 vfps::meshindex_t y0,
                                 vfps::SourceMap::hi* myhinfo)
{
    // gridpoint matrix used for interpolation
    std::unique_ptr<hi[]> ph1D(new hi[_ip]);
    std::unique_ptr<hi*[]> ph(new hi*[_it]);
    for (uint_fast8_t i=0; i<_it;i++) {
        ph[i] = &ph1D[i*_it];
    }

    // arrays of interpolation coefficients
    std::unique_ptr<interpol_t[]> icq(new interpol_t[_it]);
    std::unique_ptr<interpol_t[]> icp(new interpol_t[_it]);

    std::unique_ptr<interpol_t[]> smc(new interpol_t[_ip]);

    // new coordinates of grid point x0,y0.
    meshaxis_t x1r = (_cos_dt*_axis[0]->at(x0)-_sin_dt*_axis[1]->at(y0))
                    / _axis[0]->delta()
                    + _axis[0]->zerobin();
    meshaxis_t y1r = (_sin_dt*_axis[0]->at(x0)+_cos_dt*_axis[1]->at(y0))
                    / _axis[1]->delta()
                    + _axis[1]->zerobin();

    // new coordinates floating point part
    interpol_t xf = std::modf(x1r, &x1r);
    interpol_t yf = std::modf(y1r, &y1r);

    // new coordinates integer parts
    meshindex_t x1 = static_cast<meshindex_t>(x1r);
    meshindex_t y1 = static_cast<meshindex_t>(y1r);

    if (x1 <  _xsize && y1 < _ysize) {
        // create vectors containing interpolation coefficiants
        calcCoefficiants(icq.get(),xf,_it);
        calcCoefficiants(icp.get(),yf,_it);

        /*  Assemble interpolation
         * (using uint_fast16_t so product won't overflow)
         */
        for (uint_fast16_t smq=0; smq<_it; smq++) {
            for (uint_fast16_t smp=0; smp<_it; smp++){
                smc[smp*_it+smq] = icq[smp]*icp[smq];
            }
        }


        // renormlize to minimize rounding errors
        // renormalize(smc.size(),smc.data());

        // write source map
        for (meshindex_t j1=0; j1<_it; j1++) {
             meshindex_t j0 = y1+j1-(_it-1)/2;
            for (meshindex_t i1=0; i1<_it; i1++) {
                 meshindex_t i0 = x1+i1-(_it-1)/2;
                if(i0< _xsize && j0 < _ysize ){
                    ph[i1][j1].index = i0*_ysize+j0;
                    ph[i1][j1].weight = smc[i1*_it+j1];
                } else {
                    ph[i1][j1] = {0,0};
                }
                myhinfo[i1*_it+j1] = ph[i1][j1];
            }
        }
    } else {
        for (std::remove_const<decltype(_ip)>::type i=0; i<_ip; i++) {
            myhinfo[i] = {0,0};
        }
    }
}


#if INOVESA_USE_OPENCL == 1
void vfps::RotationMap::genCode4SM4sat()
{
    _cl_code+= R"(
    __kernel void applySM4sat(    const __global data_t* src,
                                const __global hi* sm,
                                __global data_t* dst)
    {
        data_t value = 0;
        const uint i = get_global_id(0);
        const uint offset = i*16;
    )";
    if (std::is_same<vfps::meshdata_t,float>::value) {
        _cl_code += R"(
            data_t ceil=0.0f;
            data_t flor=1.0f;
        )";
    }
    if (std::is_same<vfps::meshdata_t,double>::value) {
        _cl_code += R"(
            data_t ceil=0.0;
            data_t flor=1.0;
        )";
    }
    _cl_code += R"(
        data_t tmp;
        value += mult(src[sm[offset].src],sm[offset].weight);
        value += mult(src[sm[offset+1].src],sm[offset+1].weight);
        value += mult(src[sm[offset+2].src],sm[offset+2].weight);
        value += mult(src[sm[offset+3].src],sm[offset+3].weight);
        value += mult(src[sm[offset+4].src],sm[offset+4].weight);
        tmp = src[sm[offset+5].src];
        value += mult(tmp,sm[offset+5].weight);
        ceil = max(ceil,tmp);
        flor = min(flor,tmp);
        tmp = src[sm[offset+6].src];
        value += mult(tmp,sm[offset+6].weight);
        ceil = max(ceil,tmp);
        flor = min(flor,tmp);
        value += mult(src[sm[offset+7].src],sm[offset+7].weight);
        value += mult(src[sm[offset+8].src],sm[offset+8].weight);
        tmp = src[sm[offset+9].src];
        value += mult(tmp,sm[offset+9].weight);
        ceil = max(ceil,tmp);
        flor = min(flor,tmp);
        tmp = src[sm[offset+10].src];
        value += mult(tmp,sm[offset+10].weight);
        ceil = max(ceil,tmp);
        flor = min(flor,tmp);
        value += mult(src[sm[offset+11].src],sm[offset+11].weight);
        value += mult(src[sm[offset+12].src],sm[offset+12].weight);
        value += mult(src[sm[offset+13].src],sm[offset+13].weight);
        value += mult(src[sm[offset+14].src],sm[offset+14].weight);
        value += mult(src[sm[offset+15].src],sm[offset+15].weight);
        dst[i] = clamp(value,flor,ceil);
})";
}

void vfps::RotationMap::genCode4Rotation()
{
    _cl_code += R"(
    __kernel void applyRotation(    const __global data_t* src,
                                    const float2 rot,
                                    const int2 imgSize,
                                    const float2 zerobin,
                                    __global data_t* dst)
    {
        const int x = get_global_id(0);
        const int y = get_global_id(1);

        const data_t srcx = rot.x*(x-zerobin.x)-rot.y*(y-zerobin.y)+zerobin.x;
        const data_t srcy = rot.y*(x-zerobin.x)+rot.x*(y-zerobin.y)+zerobin.y;
        const int xi = floor(srcx);
        const int yi = floor(srcy);
        const data_t xf = srcx - xi;
        const data_t yf = srcy - yi;
    )";

    switch (_it) {
    case InterpolationType::quadratic:
        _cl_code += R"(
        const data3_t icq = (data3_t)(xf*(xf-1)/2,(1-xf*xf),xf*(xf+1)/2);
        const data3_t icp = (data3_t)(yf*(yf-1)/2,(1-yf*yf),yf*(yf+1)/2);
        )";
        break;
    case InterpolationType::cubic:
            _cl_code += R"(
            const data4_t icq = (data4_t)(
                                    (xf  )*(xf-1)*(xf-2)/(-6),
                                    (xf+1)*(xf-1)*(xf-2)/( 2),
                                    (xf+1)*(xf  )*(xf-2)/(-2),
                                    (xf+1)*(xf  )*(xf-1)/( 6)
                                );

            const data4_t icp = (data4_t)(
                                    (yf  )*(yf-1)*(yf-2)/(-6),
                                    (yf+1)*(yf-1)*(yf-2)/( 2),
                                    (yf+1)*(yf  )*(yf-2)/(-2),
                                    (yf+1)*(yf  )*(yf-1)/( 6)
                                );
            )";
            break;
    }

    if (_clamp) {
    _cl_code += R"(
        data_t hi = max(max(src[(xi  )*imgSize.y+yi],src[(xi  )*imgSize.y+yi+1]),
                       max(src[(xi+1)*imgSize.y+yi],src[(xi+1)*imgSize.y+yi+1]));
        data_t lo = min(min(src[(xi  )*imgSize.y+yi],src[(xi  )*imgSize.y+yi+1]),
                       min(src[(xi+1)*imgSize.y+yi],src[(xi+1)*imgSize.y+yi+1]));
        dst[x*imgSize.y+y] = clamp(
        )";
    } else {
        _cl_code += R"(
            dst[x*imgSize.y+y] =
            )";
    }


    switch (_it) {
    case InterpolationType::quadratic:
        _cl_code += R"(
                icq.s0*(
                    +icp.s0*src[(xi-1)*imgSize.y+yi-1]
                    +icp.s1*src[(xi-1)*imgSize.y+yi  ]
                    +icp.s2*src[(xi-1)*imgSize.y+yi+1]
                ) +
                icq.s1*(
                    +icp.s0*src[(xi  )*imgSize.y+yi-1]
                    +icp.s1*src[(xi  )*imgSize.y+yi  ]
                    +icp.s2*src[(xi  )*imgSize.y+yi+1]
                ) +
                icq.s2*(
                    +icp.s0*src[(xi+1)*imgSize.y+yi-1]
                    +icp.s1*src[(xi+1)*imgSize.y+yi  ]
                    +icp.s2*src[(xi+1)*imgSize.y+yi+1]
                )
        )";
        break;
    case InterpolationType::cubic:
        _cl_code += R"(
                icq.s0*(
                    +icp.s0*src[(xi-1)*imgSize.y+yi-1]
                    +icp.s1*src[(xi-1)*imgSize.y+yi  ]
                    +icp.s2*src[(xi-1)*imgSize.y+yi+1]
                    +icp.s3*src[(xi-1)*imgSize.y+yi+2]
                ) +
                icq.s1*(
                    +icp.s0*src[(xi  )*imgSize.y+yi-1]
                    +icp.s1*src[(xi  )*imgSize.y+yi  ]
                    +icp.s2*src[(xi  )*imgSize.y+yi+1]
                    +icp.s3*src[(xi  )*imgSize.y+yi+2]
                ) +
                icq.s2*(
                    +icp.s0*src[(xi+1)*imgSize.y+yi-1]
                    +icp.s1*src[(xi+1)*imgSize.y+yi  ]
                    +icp.s2*src[(xi+1)*imgSize.y+yi+1]
                    +icp.s3*src[(xi+1)*imgSize.y+yi+2]
                ) +
                icq.s3*(
                    +icp.s0*src[(xi+2)*imgSize.y+yi-1]
                    +icp.s1*src[(xi+2)*imgSize.y+yi  ]
                    +icp.s2*src[(xi+2)*imgSize.y+yi+1]
                    +icp.s3*src[(xi+2)*imgSize.y+yi+2]
                )
        )";
        break;
    }
    if (_clamp) {
    _cl_code += R"(
        ,lo,hi);})";
    } else {
    _cl_code += R"(
        ;})";
    }
}
#endif // INOVESA_USE_OPENCL
