/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2013-2016: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2016: Karlsruhe Institute of Technology                 *
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

#include "SM/RotationMap.hpp"

vfps::RotationMap::RotationMap(std::shared_ptr<PhaseSpace> in,
                               std::shared_ptr<PhaseSpace> out,
                               const meshindex_t xsize,
                               const meshindex_t ysize,
                               const meshaxis_t angle,
                               const InterpolationType it,
                               const bool interpol_clamped,
                               const RotationCoordinates rt,
                               const size_t rotmapsize) :
    SourceMap(in,out,xsize,ysize,size_t(rotmapsize)*it*it,it*it,it),
    _rotmapsize(rotmapsize),
    _clamp(interpol_clamped),
    _rt(rt),
    _cos_dt(cos(-angle)),
    _sin_dt(sin(-angle))
{
    if (_rotmapsize == 0) {
        #ifdef INOVESA_USE_CL
        if (OCLH::active) {
            rot = {{float(_cos_dt),float(_sin_dt)}};
            imgsize = {{cl_int(_xsize),cl_int(_ysize)}};

            genCode4Rotation();
            _cl_prog  = OCLH::prepareCLProg(_cl_code);

            applySM = cl::Kernel(_cl_prog, "applyRotation");
            applySM.setArg(0, _in->data_buf);
            applySM.setArg(1, imgsize);
            applySM.setArg(2, rot);
            applySM.setArg(3, _out->data_buf);
        } else
        #endif // INOVESA_USE_CL
        {
            if (_clamp) {
                notClampedMessage();
            }
        }
    } else {
        meshindex_t maxx;
        if (_rotmapsize == _size) {
            maxx = _xsize;
        } else {
            maxx = _xsize/2;
        }
        for (meshindex_t q_i=0; q_i< maxx; q_i++) {
            for(meshindex_t p_i=0; p_i< _ysize; p_i++) {
                genHInfo(q_i,p_i,&_hinfo[(q_i*ysize+p_i)*_ip]);
            }
        }
        #ifdef INOVESA_USE_CL
        if (OCLH::active) {
            _hi_buf = cl::Buffer(OCLH::context,
                                 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                                 sizeof(hi)*_ip*_rotmapsize,
                                 _hinfo);
            if (_clamp) {
                if (it == InterpolationType::cubic) {
                    if (_rotmapsize == _size) {
                        genCode4SM4sat();
                        _cl_prog  = OCLH::prepareCLProg(_cl_code);

                        applySM = cl::Kernel(_cl_prog, "applySM4sat");
                        applySM.setArg(0, _in->data_buf);
                        applySM.setArg(1, _hi_buf);
                        applySM.setArg(2, _out->data_buf);
                    }
                }
            } else {
                genCode4SM1D();
                _cl_prog  = OCLH::prepareCLProg(_cl_code);

                applySM = cl::Kernel(_cl_prog, "applySM1D");
                applySM.setArg(0, _in->data_buf);
                applySM.setArg(1, _hi_buf);
                applySM.setArg(2, _ip);
                applySM.setArg(3, _out->data_buf);
            }
        }
        #endif // INOVESA_USE_CL
    }
}

vfps::RotationMap::~RotationMap()
#ifdef INOVESA_ENABLE_CLPROFILING
{
    saveTimings("RotationMap");
}
#else
= default;
#endif // INOVESA_ENABLE_CLPROFILING

void vfps::RotationMap::apply()
{
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        #ifdef INOVESA_SYNC_CL
        _in->syncCLMem(clCopyDirection::cpu2dev);
        #endif // INOVESA_SYNC_CL
        if (_rotmapsize == 0) {
             // stay away from mesh borders
            OCLH::enqueueNDRangeKernel (
                        applySM,
                        cl::NDRange(1,1),
                        cl::NDRange(_xsize-_it+1,_ysize-_it+1));
        } else {
            OCLH::enqueueNDRangeKernel (
                        applySM,
                        cl::NullRange,
                        cl::NDRange(_rotmapsize));
        }
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

        if (_rotmapsize == 0) {
            for (meshindex_t q_i=0; q_i< _xsize/2; q_i++) {
                for(meshindex_t p_i=0; p_i< _ysize; p_i++) {
                    meshindex_t i = q_i*_ysize+p_i;
                    data_out[i] = 0;
                    data_out[_size-1-i] = 0;
                    genHInfo(q_i,p_i,_hinfo);
                    for (uint_fast8_t j=0; j<_ip; j++) {
                        hi h = _hinfo[j];
                        data_out[i] += data_in[h.index]*static_cast<meshdata_t>(h.weight);
                        data_out[_size-1-i] += data_in[_size-1-h.index]*static_cast<meshdata_t>(h.weight);
                    }
                }
            }
        } else {
            for (meshindex_t i=0; i< _rotmapsize; i++) {
                data_out[i] = 0;
                if (_rotmapsize == _size/2) {
                    data_out[_size-1-i] = 0;
                }
                for (uint_fast8_t j=0; j<_ip; j++) {
                    hi h = _hinfo[i*_ip+j];
                    data_out[i] += data_in[h.index]*static_cast<meshdata_t>(h.weight);
                    if (_rotmapsize == _size/2) {
                        data_out[_size-1-i] += data_in[_size-1-h.index]*static_cast<meshdata_t>(h.weight);
                    }
                }
                if (_clamp) {
                    // handle overshooting
                    meshdata_t ceil=std::numeric_limits<meshdata_t>::min();
                    meshdata_t flor=std::numeric_limits<meshdata_t>::max();
                    for (size_t x=1; x<=2; x++) {
                        for (size_t y=1; y<=2; y++) {
                            ceil = std::max(ceil,data_in[_hinfo[i*_ip+x*_it+y].index]);
                            flor = std::min(flor,data_in[_hinfo[i*_ip+x*_it+y].index]);
                        }
                    }
                    data_out[i] = std::max(std::min(ceil,data_out[i]),flor);

                    if (_rotmapsize == _size/2) {
                        ceil=std::numeric_limits<meshdata_t>::min();
                        flor=std::numeric_limits<meshdata_t>::max();
                        for (size_t x=1; x<=2; x++) {
                            for (size_t y=1; y<=2; y++) {
                                ceil = std::max(ceil,data_in[_size-1-_hinfo[i*_ip+x*_it+y].index]);
                                flor = std::min(flor,data_in[_size-1-_hinfo[i*_ip+x*_it+y].index]);
                            }
                        }
                        data_out[_size-1-i] = std::max(std::min(ceil,data_out[_size-1-i]),flor);
                    }
                }
            }
        }
    }
}

vfps::PhaseSpace::Position
vfps::RotationMap::apply(const PhaseSpace::Position pos) const
{
    PhaseSpace::Position rv;
    rv.x = _cos_dt*meshaxis_t(pos.x-(_xsize-1)/2.0)
         + _sin_dt*meshaxis_t(pos.y-(_ysize-1)/2.0)
         + meshaxis_t((_xsize-1)/2.0);
    rv.y = _cos_dt*meshaxis_t(pos.y-(_ysize-1)/2.0)
         - _sin_dt*meshaxis_t(pos.x-(_xsize-1)/2.0)
         + meshaxis_t((_ysize-1)/2.0);
    return rv;
}

void vfps::RotationMap::genHInfo(vfps::meshindex_t q_i,
                                 vfps::meshindex_t p_i,
                                 vfps::SourceMap::hi* myhinfo)
{
    // gridpoint matrix used for interpolation
    hi* ph1D = new hi[_ip];
    hi** ph = new hi*[_it];
    for (uint_fast8_t i=0; i<_it;i++) {
        ph[i] = &ph1D[i*_it];
    }

    // arrays of interpolation coefficients
    interpol_t* icq = new interpol_t[_it];
    interpol_t* icp = new interpol_t[_it];

    interpol_t* smc = new interpol_t[_ip];

    // Cell of inverse image (qp,pp) of grid point i,j.
    meshaxis_t qp; //q', backward mapping
    meshaxis_t pp; //p'
    // interpolation type specific q and p coordinates
    meshaxis_t pcoord;
    meshaxis_t qcoord;
    meshaxis_t qq_int;
    meshaxis_t qp_int;
    //Scaled arguments of interpolation functions:
    meshindex_t xi; //meshpoint smaller q'
    meshindex_t yi; //numper of lower mesh point from p'
    interpol_t xf; //distance from id
    interpol_t yf; //distance of p' from lower mesh point
    switch (_rt) {
    case RotationCoordinates::mesh:
        qp = _cos_dt*meshaxis_t(q_i-(_xsize-1)/2.0)
                - _sin_dt*meshaxis_t(p_i-(_ysize-1)/2.0)
                + meshaxis_t((_xsize-1)/2.0);
        pp = _sin_dt*meshaxis_t(q_i-(_xsize-1)/2.0)
                + _cos_dt*meshaxis_t(p_i-(_ysize-1)/2.0)
                +meshaxis_t((_ysize-1)/2.0);
        qcoord = qp;
        pcoord = pp;
        break;
    case RotationCoordinates::norm_0_1:
        qp = _cos_dt*meshaxis_t((q_i-(_xsize-1)/2.0)/(_xsize-1))
           - _sin_dt*meshaxis_t((p_i-(_ysize-1)/2.0)/(_ysize-1));
        pp = _sin_dt*meshaxis_t((q_i-(_xsize-1)/2.0)/(_xsize-1))
           + _cos_dt*meshaxis_t((p_i-(_ysize-1)/2.0)/(_ysize-1));
        qcoord = (qp+meshaxis_t(0.5))*meshaxis_t(_xsize-1);
        pcoord = (pp+meshaxis_t(0.5))*meshaxis_t(_ysize-1);
        break;
    case RotationCoordinates::norm_pm1:
    default:
        qp = _cos_dt*meshaxis_t(2*int(q_i)-int(_xsize-1))
                    /meshaxis_t(_xsize-1)
           - _sin_dt*meshaxis_t(2*int(p_i)-int(_ysize-1))
                    /meshaxis_t(_ysize-1);

        pp = _sin_dt*meshaxis_t(2*int(q_i)-int(_xsize-1))
                    /meshaxis_t(_xsize-1)
           + _cos_dt*meshaxis_t(2*int(p_i)-int(_ysize-1))
                    /meshaxis_t(_ysize-1);
        qcoord = (qp+meshaxis_t(1))*meshaxis_t(_xsize-1)/meshaxis_t(2);
        pcoord = (pp+meshaxis_t(1))*meshaxis_t(_ysize-1)/meshaxis_t(2);
        break;
    }
    xf = modf(qcoord, &qq_int);
    yf = modf(pcoord, &qp_int);
    xi = qq_int;
    yi = qp_int;

    if (xi <  _xsize && yi < _ysize) {
        // create vectors containing interpolation coefficiants
        calcCoefficiants(icq,xf,_it);
        calcCoefficiants(icp,yf,_it);

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

        // write heritage map
        for (meshindex_t j1=0; j1<_it; j1++) {
             meshindex_t j0 = yi+j1-(_it-1)/2;
            for (meshindex_t i1=0; i1<_it; i1++) {
                 meshindex_t i0 = xi+i1-(_it-1)/2;
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
        for (uint_fast8_t i=0; i<_ip; i++) {
            myhinfo[i] = {0,0};
        }
    }

    delete [] icp;
    delete [] icq;
    delete [] smc;

    delete [] ph;
    delete [] ph1D;
}


#ifdef INOVESA_USE_CL
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
    #if FXP_FRACPART < 31
    if (std::is_same<vfps::meshdata_t,vfps::fixp32>::value) {
        _cl_code += R"(
            data_t ceil=INT_MIN;
            data_t flor=INT_MAX;
        )";
    }
    #endif
    if (std::is_same<vfps::meshdata_t,vfps::fixp64>::value) {
        _cl_code += R"(
            data_t ceil=LONG_MIN;
            data_t flor=LONG_MAX;
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
                                    const int2 imgSize,
                                    const float2 rot,
                                    __global data_t* dst)
    {
        const int x = get_global_id(0)+1;
        const int y = get_global_id(1)+1;

        const data_t srcx = rot.x*(x-(imgSize.x+1)/2)-rot.y*(y-(imgSize.y+1)/2)
                                +(imgSize.x+1)/2;
        const data_t srcy = rot.y*(x-(imgSize.x+1)/2)+rot.x*(y-(imgSize.y+1)/2)
                                +(imgSize.y+1)/2;
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
#endif // INOVESA_USE_CL
