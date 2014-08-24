#include "CL/CLProgApplyHM.hpp"

void prepareCLProgApplyHM()
{
	const char* code = R"(
		typedef struct {
			uint src;
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
		)";

	cl::Program::Sources source(1,std::make_pair(code,strlen(code)));
	CLProgApplyHM::p = cl::Program(OCLH::context, source);
	CLProgApplyHM::p.build(OCLH::devices);
#ifdef DEBUG
	std::cout	<< "Build Log:\t "
				<< CLProgApplyHM::p.getBuildInfo<CL_PROGRAM_BUILD_LOG>(OCLH::devices[0])
				<< std::endl;
#endif // DEBUG
}

cl::Program CLProgApplyHM::p;
