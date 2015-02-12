#include "HM/HeritageMap.hpp"

vfps::HeritageMap::HeritageMap(meshdata_t* in, meshdata_t* out,
							   size_t xsize, size_t ysize) :
	_data1D(in),
	_data1D_rotated(out),
	_heritage_map(new std::array<hi,it*it>*[xsize]),
	_heritage_map1D(new std::array<hi,it*it>[xsize*ysize]()),
	_size(xsize*ysize),
	_xsize(xsize),
	_ysize(ysize)
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
#ifdef FR_USE_CL
	OCLH::queue.enqueueNDRangeKernel (
				applyHM2D,
				cl::NullRange,
				cl::NDRange(size(0),size(1)));
#ifdef CL_VERSION_1_2
	OCLH::queue.enqueueBarrierWithWaitList();
#else // CL_VERSION_1_2
	OCLH::queue.enqueueBarrier();
#endif // CL_VERSION_1_2
	cl::size_t<3> null3d;
	cl::size_t<3> imgsize;
	imgsize[0] = size(0);
	imgsize[1] = size(1);
	imgsize[2] = 1;
	OCLH::queue.enqueueCopyImage(_data_rotated_buf,
								 _data_buf,
								 null3d,null3d,imgsize);
#else // FR_USE_CL
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
#endif // FR_USE_CL
}
