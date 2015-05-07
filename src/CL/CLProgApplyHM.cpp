/******************************************************************************/
/* Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   */
/* Copyright (c) 2014-2015: Patrik Schönfeldt                                 */
/*                                                                            */
/* This file is part of Inovesa.                                              */
/* Inovesa is free software: you can redistribute it and/or modify            */
/* it under the terms of the GNU General Public License as published by       */
/* the Free Software Foundation, either version 3 of the License, or          */
/* (at your option) any later version.                                        */
/*                                                                            */
/* Inovesa is distributed in the hope that it will be useful,                 */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           */
/******************************************************************************/

#include "CL/CLProgApplyHM.hpp"

void prepareCLProgApplyHM()
{
	std::string code;

	if (std::is_same<vfps::meshdata_t,float>::value ||
		std::is_same<vfps::meshdata_t,double>::value) {
		code = R"(
		typedef struct {
			uint src;
			float weight;
		} hi;

		__kernel void applyHM1D(const __global float* src,
								const __global hi* hm,
								const uint hm_len,
								__global float* dst)
		{
			float value = 0;
			const uint i = get_global_id(0);
			const uint offset = i*hm_len;
			for (uint j=0; j<hm_len; j++)
			{
				value += src[hm[offset+j].src] * hm[offset+j].weight;
			}
			dst[i] = value;
		}

		__kernel void applyHM4sat(	const __global float* src,
									const __global hi* hm,
									__global float* dst)
		{
			float value = 0;
			const uint i = get_global_id(0);
			const uint offset = i*16;
			float ceil=0;
			float flor=1;
			float tmp;
			value += src[hm[offset].src] * hm[offset].weight;
			value += src[hm[offset+1].src] * hm[offset+1].weight;
			value += src[hm[offset+2].src] * hm[offset+2].weight;
			value += src[hm[offset+3].src] * hm[offset+3].weight;
			value += src[hm[offset+4].src] * hm[offset+4].weight;
			tmp = src[hm[offset+5].src];
			value += tmp * hm[offset+5].weight;
			ceil = max(ceil,tmp);
			flor = min(flor,tmp);
			tmp = src[hm[offset+6].src];
			value += tmp * hm[offset+6].weight;
			ceil = max(ceil,tmp);
			flor = min(flor,tmp);
			value += src[hm[offset+7].src] * hm[offset+7].weight;
			value += src[hm[offset+8].src] * hm[offset+8].weight;
			tmp = src[hm[offset+9].src];
			value += tmp * hm[offset+9].weight;
			ceil = max(ceil,tmp);
			flor = min(flor,tmp);
			tmp = src[hm[offset+10].src];
			value += tmp * hm[offset+10].weight;
			ceil = max(ceil,tmp);
			flor = min(flor,tmp);
			value += src[hm[offset+11].src] * hm[offset+11].weight;
			value += src[hm[offset+12].src] * hm[offset+12].weight;
			value += src[hm[offset+13].src] * hm[offset+13].weight;
			value += src[hm[offset+14].src] * hm[offset+14].weight;
			value += src[hm[offset+15].src] * hm[offset+15].weight;

			dst[i] = max(min(ceil,value),flor);
		}
		)";
	} else {
		std::stringstream fxp_fracpart;
		fxp_fracpart << FXP_FRACPART;
		if (std::is_same<vfps::meshdata_t,vfps::fixp64>::value) {
			code = R"(
			typedef long meshdata_t;
			typedef long interpol_t;
			typedef long integral_t;
			)";
		} else {
			code = R"(
			typedef int meshdata_t;
			typedef int interpol_t;
			typedef long integral_t;
			)";
		}

		code += R"(
		typedef struct {
			uint src;
			interpol_t weight;
		} hi;

		__kernel void applyHM1D(const __global meshdata_t* src,
								const __global hi* hm,
								const uint hm_len,
								__global meshdata_t* dst)
		{
			integral_t value = 0;
			const uint i = get_global_id(0);
			const uint offset = i*hm_len;
			for (uint j=0; j<hm_len; j++)
			{
				value += (integral_t)(src[hm[offset+j].src])
						* (integral_t)(hm[offset+j].weight);
			}
			dst[i] = value >>
		)";
		code += fxp_fracpart.str();
		code += R"(
					;
		}
		__kernel void applyHM4sat(	const __global meshdata_t* src,
									const __global hi* hm,
									__global meshdata_t* dst)
		{
			integral_t value = 0;
			meshdata_t result;
			const uint i = get_global_id(0);
			const uint offset = i*16;
		)";
		if (std::is_same<vfps::meshdata_t,vfps::fixp64>::value) {
			code += R"(
				meshdata_t ceil=LONG_MIN;
				meshdata_t flor=LONG_MAX;
			)";
		} else {
			code += R"(
				meshdata_t ceil=INT_MIN;
				meshdata_t flor=INT_MAX;
			)";
		}
		code += R"(
			meshdata_t tmp;
			value += (integral_t)(src[hm[offset].src])
					* (integral_t)(hm[offset].weight);
			value += (integral_t)(src[hm[offset+1].src])
					* (integral_t)(hm[offset+1].weight);
			value += (integral_t)(src[hm[offset+2].src])
					* (integral_t)(hm[offset+2].weight);
			value += (integral_t)(src[hm[offset+3].src])
					* (integral_t)(hm[offset+3].weight);
			value += (integral_t)(src[hm[offset+4].src])
					* (integral_t)(hm[offset+4].weight);
			tmp = src[hm[offset+5].src];
			value += (integral_t)(tmp)
					* (integral_t)(hm[offset+5].weight);
			ceil = max(ceil,tmp);
			flor = min(flor,tmp);
			tmp = src[hm[offset+6].src];
			value += (integral_t)(tmp)
					* (integral_t)(hm[offset+6].weight);
			ceil = max(ceil,tmp);
			flor = min(flor,tmp);
			value += (integral_t)(src[hm[offset+7].src])
					* (integral_t)(hm[offset+7].weight);
			value += (integral_t)(src[hm[offset+8].src])
					* (integral_t)(hm[offset+8].weight);
			tmp = src[hm[offset+9].src];
			value += (integral_t)(tmp)
					* (integral_t)(hm[offset+9].weight);
			ceil = max(ceil,tmp);
			flor = min(flor,tmp);
			tmp = src[hm[offset+10].src];
			value += (integral_t)(tmp)
					* (integral_t)(hm[offset+10].weight);
			ceil = max(ceil,tmp);
			flor = min(flor,tmp);
			value += (integral_t)(src[hm[offset+11].src])
					* (integral_t)(hm[offset+11].weight);
			value += (integral_t)(src[hm[offset+12].src])
					* (integral_t)(hm[offset+12].weight);
			value += (integral_t)(src[hm[offset+13].src])
					* (integral_t)(hm[offset+13].weight);
			value += (integral_t)(src[hm[offset+14].src])
					* (integral_t)(hm[offset+14].weight);
			value += (integral_t)(src[hm[offset+15].src])
					* (integral_t)(hm[offset+15].weight);

			result = value >>
		)";
		code += fxp_fracpart.str();
		code += R"(
					;
			dst[i] = max(min(ceil,result),flor);
		}
		)";
	}

	cl::Program::Sources source(1,std::make_pair(code.c_str(),code.length()));
	CLProgApplyHM::p = cl::Program(OCLH::context, source);
    try {
        CLProgApplyHM::p.build(OCLH::devices);
    } catch (cl::Error &e) {
        e.what();
    #ifdef DEBUG
        std::cout	<< "Build Log:\t "
                    << CLProgApplyHM::p.getBuildInfo<CL_PROGRAM_BUILD_LOG>(OCLH::devices[0])
                    << std::endl;
    #endif // DEBUG
    }
}

cl::Program CLProgApplyHM::p;
