#include "CL/CLProgApplyHM.hpp"

void prepareCLProgApplyHM()
{
	const char* code = R"(
		typedef struct {
			uint src;
			int2 src2d;
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

		__kernel void applyHM2D(__read_only image2d_t src,
								const __global hi* hm,
								uint img_height,
								const uint hm_len,
								__write_only image2d_t dst)
		{
			float4 value = 0;
			const uint x = get_global_id(0);
			const uint y = get_global_id(1);
			int2 coords = (int2)(x,y);
			const sampler_t sampler
					= CLK_NORMALIZED_COORDS_FALSE
					| CLK_ADDRESS_NONE
					| CLK_FILTER_NEAREST;
			const uint i = x*img_height+y;
			const uint offset = i*hm_len;
			for (uint j=0; j<hm_len; j++)
			{
				value += read_imagef(src,sampler,hm[offset+j].src2d)
						* hm[offset+j].weight;
			}
			write_imagef(dst,coords,value);
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

	cl::Program::Sources source(1,std::make_pair(code,strlen(code)));
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
