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

	std::stringstream fxp_fracpart;
	fxp_fracpart << FXP_FRACPART;
	std::stringstream fxp_intmask;
	fxp_intmask << std::hex <<~vfps::fixp64(-1.0).internalValue();

	code =	"\t__constant int fracpart="+fxp_fracpart.str()+";\n"
			"\t__constant long intmask=0x"+fxp_intmask.str()+";\n";
	if (std::is_same<vfps::meshdata_t,float>::value) {
	code += R"(
	typedef float data_t;
	float mult(float x, float y) { return x*y; }
	)";
	}
	if (std::is_same<vfps::meshdata_t,double>::value) {
	code += R"(
	typedef double data_t;
	double mult(double x, double y) { return x*y; }
	)";
	}
	#if FXP_FRACPART < 31
	if (std::is_same<vfps::meshdata_t,vfps::fixp32>::value) {
	code += R"(
	typedef int data_t;
	int mult(int x, int y) {
		return ((long)(x)*(long)(y))>>fracpart;
	}
	)";
	}
	#endif
	if (std::is_same<vfps::meshdata_t,vfps::fixp64>::value) {
	code += R"(
	typedef long data_t;
	ulong umult(ulong x, ulong y) {
		ulong a = x >> fracpart;
		ulong b = (x & intmask);
		ulong c = y >> fracpart;
		ulong d = (y & intmask);
		return ((d*b) >> fracpart) + (d*a) + (c*b) + ((c*a) << fracpart);
	}
	long mult(long x,long y) {
		long val = umult(abs(x),abs(y));
		if ((x<0 && y>0)||(x>0 && y<0))
			return -val;
		return val;
	}

	)";
	}

	code += R"(
	typedef struct {
		uint src;
		data_t weight;
	} hi;

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
/*
	code += R"(
	__kernel void applyHM4sat(	const __global data_t* src,
								const __global hi* hm,
								__global data_t* dst)
	{
		data_t value = 0;
		data_t result;
		const uint i = get_global_id(0);
		const uint offset = i*16;
	)";
	if (std::is_same<vfps::meshdata_t,float>::value) {
		code += R"(
			data_t ceil=0.0f;
			data_t flor=1.0f;
		)";
	}
	if (std::is_same<vfps::meshdata_t,double>::value) {
		code += R"(
			data_t ceil=0.0;
			data_t flor=1.0;
		)";
	}
	#if FXP_FRACPART < 31
	if (std::is_same<vfps::meshdata_t,vfps::fixp32>::value) {
		code += R"(
			data_t ceil=INT_MIN;
			data_t flor=INT_MAX;
		)";
	}
	#endif
	if (std::is_same<vfps::meshdata_t,vfps::fixp64>::value) {
		code += R"(
			data_t ceil=LONG_MIN;
			data_t flor=LONG_MAX;
		)";
	}
	code += R"(
		 data_t tmp;
		value += (data_t)(src[hm[offset].src])
				* (data_t)(hm[offset].weight);
		value += (data_t)(src[hm[offset+1].src])
				* (data_t)(hm[offset+1].weight);
		value += (data_t)(src[hm[offset+2].src])
				* (data_t)(hm[offset+2].weight);
		value += (data_t)(src[hm[offset+3].src])
				* (data_t)(hm[offset+3].weight);
		value += (data_t)(src[hm[offset+4].src])
				* (data_t)(hm[offset+4].weight);
		tmp = src[hm[offset+5].src];
		value += (data_t)(tmp)
				* (data_t)(hm[offset+5].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		tmp = src[hm[offset+6].src];
		value += (data_t)(tmp)
				* (data_t)(hm[offset+6].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		value += ( data_t)(src[hm[offset+7].src])
				* ( data_t)(hm[offset+7].weight);
		value += ( data_t)(src[hm[offset+8].src])
				* ( data_t)(hm[offset+8].weight);
		tmp = src[hm[offset+9].src];
		value += ( data_t)(tmp)
				* ( data_t)(hm[offset+9].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		tmp = src[hm[offset+10].src];
		value += ( data_t)(tmp)
				* ( data_t)(hm[offset+10].weight);
		ceil = max(ceil,tmp);
		flor = min(flor,tmp);
		value += ( data_t)(src[hm[offset+11].src])
				* ( data_t)(hm[offset+11].weight);
		value += ( data_t)(src[hm[offset+12].src])
				* ( data_t)(hm[offset+12].weight);
		value += ( data_t)(src[hm[offset+13].src])
				* ( data_t)(hm[offset+13].weight);
		value += ( data_t)(src[hm[offset+14].src])
				* ( data_t)(hm[offset+14].weight);
		value += ( data_t)(src[hm[offset+15].src])
				* ( data_t)(hm[offset+15].weight);
		)";

	if (std::is_same<vfps::meshdata_t,vfps::fixp64>::value
		#if FXP_FRACPART < 31
		|| std::is_same<vfps::meshdata_t,vfps::fixp32>::value
		#endif
		) {
		code +="result = value >> " + fxp_fracpart.str() + ';';
	} else {
		code +="result = value;";
	}
	code += R"(
		dst[i] = max(min(ceil,result),flor);
	}
	)";
*/

	cl::Program::Sources source(1,std::make_pair(code.c_str(),code.length()));
	CLProgApplyHM::p = cl::Program(OCLH::context, source);
    try {
        CLProgApplyHM::p.build(OCLH::devices);
    } catch (cl::Error &e) {
		e.what();
	#ifdef DEBUG
	// in debug builds, CL code and build log should always be displayed
    }
	#endif  // DEBUG
		std::cout	<< "===== OpenCL Code =====\n"
					<< code << std::endl;
		std::cout	<< "===== OpenCL Build Log =====\n"
					<< CLProgApplyHM::p.getBuildInfo<CL_PROGRAM_BUILD_LOG>(OCLH::devices[0])
					<< std::endl;
	#ifndef DEBUG
	// in release builds show CL code and build log only when error occured
	}
	#endif // DEBUG
}

cl::Program CLProgApplyHM::p;
