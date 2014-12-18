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
			float4 value;
            int2 center = (int2)(imgSize.x/2,imgSize.y/2);

			int2 dstpos = (int2)(x,y);

			x -= center.x;
            y -= center.y;

			float srcx = (rot.x*(x) - rot.y*(y) + center.x);
			float srcy = (rot.y*(x) + rot.x*(y) + center.y);
			int2 srcpos = (int2)(srcx,srcy);
			float xf = srcx-srcpos.x;
			float yf = srcy-srcpos.y;

            const sampler_t sampler
                    = CLK_NORMALIZED_COORDS_FALSE
                    | CLK_ADDRESS_CLAMP_TO_EDGE
                    | CLK_FILTER_NEAREST;

			value = (
					xf*(xf-1)*(
                        +yf*(yf-1)*read_imagef(src,sampler,srcpos+(int2)(-1,-1))
                        +2*(1-yf*yf)*read_imagef(src,sampler,srcpos+(int2)(-1,0))
                        +yf*(yf+1)*read_imagef(src,sampler,srcpos+(int2)(-1,1))
					) +
					2*(1-xf*xf)*(
                        +yf*(yf-1)*read_imagef(src,sampler,srcpos+(int2)(0,-1))
                        +2*(1-yf*yf)*read_imagef(src,sampler,srcpos)
                        +yf*(yf+1)*read_imagef(src,sampler,srcpos+(int2)(0,1))
					) +
					xf*(xf+1)*(
                        +yf*(yf-1)*read_imagef(src,sampler,srcpos+(int2)(1,-1))
                        +2*(1-yf*yf)*read_imagef(src,sampler,srcpos+(int2)(1,0))
                        +yf*(yf+1)*read_imagef(src,sampler,srcpos+(int2)(1,1))
					))/4;
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
