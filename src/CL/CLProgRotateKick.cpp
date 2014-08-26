#include "CL/CLProgRotateKick.hpp"

void prepareCLProgRotateKick()
{
	const char* code = R"(
		__kernel void rotateKick (
				__write_only image2d_t dst,
				__read_only image2d_t src,
				uint2 imgSize,
				float2 rot)
		{
			int x = get_global_id(0);
			int y = get_global_id(1);
			float4 value = {0.f,0.f,0.f,0.f};
			uint2 center = imgSize/2;

			int2 dstpos = (int2)(x,y);

			if (dstpos.x == center.x && dstpos.y == center.y) value.x = 1.0f;

			x -= center.x;
			y -= center.y;

			float srcx = (rot.x*(x) - rot.y*(y) + center.x);
			float srcy = (rot.y*(x) + rot.x*(y) + center.y);
			float2 srcpos = (float2)(srcx,srcy);
			const sampler_t sampler
					= CLK_NORMALIZED_COORDS_FALSE
					| CLK_ADDRESS_CLAMP
					| CLK_FILTER_LINEAR;
			value += read_imagef(src,sampler,srcpos);
			write_imagef(dst,dstpos,value);
		}
		)";

	cl::Program::Sources source(1,std::make_pair(code,strlen(code)));
	CLProgRotateKick::p = cl::Program(OCLH::context, source);
	try {
		CLProgRotateKick::p.build(OCLH::devices);
	} catch (cl::Error &e) {
		std::cout	<< e.what() << std::endl;
		std::cout	<< "Build Log:\t "
					<< CLProgRotateKick::p.getBuildInfo<CL_PROGRAM_BUILD_LOG>(OCLH::devices[0])
					<< std::endl;
	}
}

cl::Program CLProgRotateKick::p;
