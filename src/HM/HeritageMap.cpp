#include "HM/HeritageMap.hpp"

vfps::HeritageMap::HeritageMap(PhaseSpace* in, PhaseSpace* out,
							   size_t xsize, size_t ysize) :
	_heritage_map(new std::array<hi,it*it>*[xsize]),
	_heritage_map1D(new std::array<hi,it*it>[xsize*ysize]()),
	_size(xsize*ysize),
	_xsize(xsize),
	_ysize(ysize),
	_in(in),
	_out(out)
{
	for (unsigned int i=0; i<xsize; i++) {
		_heritage_map[i] = &(_heritage_map1D[i*ysize]);
	}
}

vfps::HeritageMap::~HeritageMap()
{
	delete [] _heritage_map1D;
}

void vfps::HeritageMap::apply()
{
	#ifdef INOVESA_USE_CL
	OCLH::queue.enqueueWriteBuffer
				(_in->data_buf, CL_TRUE,
				 0,sizeof(float)*ps_xsize*ps_ysize,
				_in->getData());
	OCLH::queue.enqueueNDRangeKernel (
				applyHM,
				cl::NullRange,
				cl::NDRange(_size));
	#ifdef CL_VERSION_1_2
	OCLH::queue.enqueueBarrierWithWaitList();
	#else // CL_VERSION_1_2
	OCLH::queue.enqueueBarrier();
	#endif // CL_VERSION_1_2
	OCLH::queue.enqueueReadBuffer
				(_out->data_buf, CL_TRUE,
				 0,sizeof(float)*_size,
				_out->getData());
	#else // INOVESA_USE_CL
	for (unsigned int i=0; i< _size; i++) {
		_data1D_rotated[i] = 0;
		for (hi h: _heritage_map1D[i]) {
			_data1D_rotated[i] += _data1D[h.index]*static_cast<meshdata_t>(h.weight);
		}
		// handle overshooting
		meshdata_t ceil=std::numeric_limits<fixp32>::min();
		meshdata_t flor=std::numeric_limits<fixp32>::max();
		for (size_t x=1; x<=2; x++) {
			for (size_t y=1; y<=2; y++) {
				ceil = std::max(ceil,_data1D[_heritage_map1D[i][x*it+y].index]);
				flor = std::min(flor,_data1D[_heritage_map1D[i][x*it+y].index]);
			}
		}
		_data1D_rotated[i] = std::max(std::min(ceil,_data1D_rotated[i]),flor);
	}
	std::copy_n(_data1D_rotated,_size,_data1D);
	#endif // INOVESA_USE_CL
}

#ifdef INOVESA_USE_CL
void vfps::HeritageMap::__initOpenCL()
{
	_heritage_map1D_buf = cl::Buffer(OCLH::context,
									 CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
									 sizeof(hi)*it*it*_size,
									 _heritage_map1D);
	if (it == 4) {
		applyHM = cl::Kernel(CLProgApplyHM::p, "applyHM4sat");
		applyHM.setArg(0, _in->data_buf);
		applyHM.setArg(1, _heritage_map1D_buf);
		applyHM.setArg(2, _out->data_buf);
	} else {
		applyHM = cl::Kernel(CLProgApplyHM::p, "applyHM1D");
		applyHM.setArg(0, _in->data_buf);
		applyHM.setArg(1, _heritage_map1D_buf);
		applyHM.setArg(2, it*it);
		applyHM.setArg(3, _out->data_buf);
	}
}
#endif
