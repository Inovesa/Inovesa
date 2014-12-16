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
								__global float* dst)
		{
			const uint hm_len = 16;
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
								__write_only image2d_t dst)
		{
			const uint hm_len = 16;
			float4 value = 0;
			const uint x = get_global_id(0);
			const uint y = get_global_id(1);
			int2 coords = (int2)(x,y);
			const sampler_t sampler
					= CLK_NORMALIZED_COORDS_FALSE
					| CLK_ADDRESS_CLAMP_TO_EDGE
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
