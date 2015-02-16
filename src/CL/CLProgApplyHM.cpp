#include "CL/CLProgApplyHM.hpp"

void prepareCLProgApplyHM()
{
	const char* code;

	if (std::is_same<vfps::meshdata_t,float>::value) {
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
		code = R"(
		typedef struct {
			uint src;
			int weight;
		} hi;

		__kernel void applyHM1D(const __global int* src,
								const __global hi* hm,
								const uint hm_len,
								__global int* dst)
		{
			long value = 0;
			const uint i = get_global_id(0);
			const uint offset = i*hm_len;
			for (uint j=0; j<hm_len; j++)
			{
				value += (long)(src[hm[offset+j].src])
						* (long)(hm[offset+j].weight);
			}
			// hardcoded for 1:2:29 bit fixed point data
			dst[i] = value >> 29;
		}

		__kernel void applyHM4sat(	const __global int* src,
									const __global hi* hm,
									__global int* dst)
		{
			long value = 0;
			int result;
			const uint i = get_global_id(0);
			const uint offset = i*16;
			int ceil=INT_MIN;
			int flor=INT_MAX;
			int tmp;
			value += (long)(src[hm[offset].src])
					* (long)(hm[offset].weight);
			value += (long)(src[hm[offset+1].src])
					* (long)(hm[offset+1].weight);
			value += (long)(src[hm[offset+2].src])
					* (long)(hm[offset+2].weight);
			value += (long)(src[hm[offset+3].src])
					* (long)(hm[offset+3].weight);
			value += (long)(src[hm[offset+4].src])
					* (long)(hm[offset+4].weight);
			tmp = src[hm[offset+5].src];
			value += (long)(tmp)
					* (long)(hm[offset+5].weight);
			ceil = max(ceil,tmp);
			flor = min(flor,tmp);
			tmp = src[hm[offset+6].src];
			value += (long)(tmp)
					* (long)(hm[offset+6].weight);
			ceil = max(ceil,tmp);
			flor = min(flor,tmp);
			value += (long)(src[hm[offset+7].src])
					* (long)(hm[offset+7].weight);
			value += (long)(src[hm[offset+8].src])
					* (long)(hm[offset+8].weight);
			tmp = src[hm[offset+9].src];
			value += (long)(tmp)
					* (long)(hm[offset+9].weight);
			ceil = max(ceil,tmp);
			flor = min(flor,tmp);
			tmp = src[hm[offset+10].src];
			value += (long)(tmp)
					* (long)(hm[offset+10].weight);
			ceil = max(ceil,tmp);
			flor = min(flor,tmp);
			value += (long)(src[hm[offset+11].src])
					* (long)(hm[offset+11].weight);
			value += (long)(src[hm[offset+12].src])
					* (long)(hm[offset+12].weight);
			value += (long)(src[hm[offset+13].src])
					* (long)(hm[offset+13].weight);
			value += (long)(src[hm[offset+14].src])
					* (long)(hm[offset+14].weight);
			value += (long)(src[hm[offset+15].src])
					* (long)(hm[offset+15].weight);

			// hardcoded for 1:2:29 bit fixed point data
			result = value >> 29;

			dst[i] = max(min(ceil,result),flor);
		}
		)";
	}

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
