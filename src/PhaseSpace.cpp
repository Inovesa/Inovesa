#include "PhaseSpace.hpp"

vfps::PhaseSpace::PhaseSpace(std::array<Ruler<meshaxis_t>,2> axis) :
	_axis(axis),
	_data1D(new meshdata_t[size(0)*size(1)]()),
	_data1D_rotated(new meshdata_t[size(0)*size(1)]())
{
	_data = new meshdata_t*[size(0)];
	_data_rotated = new meshdata_t*[size(0)];
	for (unsigned int i=0; i<size(0); i++) {
		_data[i] = &(_data1D[i*size(1)]);
		_data_rotated[i] = &(_data1D_rotated[i*size(1)]);
	}
	_projection[0] = new meshdata_t[size(0)];
	_projection[1] = new meshdata_t[size(1)];
	_ws[0] = new meshdata_t[size(0)];
	_ws[1] = new meshdata_t[size(1)];

	const meshdata_t ca = 3.;
	meshdata_t dc = 1;

	const meshdata_t h03 = getDelta(0)/3.;
	_ws[0][0] = h03;
	for (unsigned int x=1; x< size(0)-1; x++){
		_ws[0][x] = h03 * (ca+dc);
		dc = -dc;
	}
	_ws[0][size(0)-1] = h03;


	const meshdata_t h13 = getDelta(1)/3.;
	_ws[1][0] = h13;
	for (unsigned int x=1; x< size(1)-1; x++){
		_ws[1][x] = h13 * (ca+dc);
		dc = -dc;
	}
	_ws[1][size(0)-1] = h13;

	img_size = {{size(0),size(1)}};
}


vfps::PhaseSpace::PhaseSpace(Ruler<meshaxis_t> axis1, Ruler<meshaxis_t> axis2) :
	PhaseSpace(std::array<Ruler<meshaxis_t>,2>{{axis1,axis2}})
{}

vfps::PhaseSpace::PhaseSpace(unsigned int xdim, meshaxis_t xmin, meshaxis_t xmax,
							 unsigned int ydim, meshaxis_t ymin, meshaxis_t ymax) :
	PhaseSpace(Ruler<meshaxis_t>(xdim,xmin,xmax),
			   Ruler<meshaxis_t>(ydim,ymin,ymax))
{}

vfps::PhaseSpace::PhaseSpace(const vfps::PhaseSpace& other) :
	PhaseSpace(other._axis)
{ std::copy_n(other._data1D,size(0)*size(1),_data1D); }

vfps::PhaseSpace::~PhaseSpace()
{
	delete [] _data;
	delete [] _data_rotated;
	delete [] _data1D;
	delete [] _data1D_rotated;
	delete [] _projection[0];
	delete [] _projection[1];
	delete [] _ws[0];
	delete [] _ws[1];
}

vfps::meshdata_t vfps::PhaseSpace::average(const unsigned int axis)
{
	if (axis == 0) {
		projectionToX();
	} else {
		projectionToY();
	}

	if (_moment[axis].size() == 0) {
		_moment[axis].resize(1);
	}

	meshdata_t avg = 0;
	for (unsigned int i=0; i<size(axis); i++) {
		avg += meshdata_t(i)*_projection[axis][i];
	}
	avg = avg/meshdata_t(size(axis));

	_moment[axis][0] = avg;

	return x(axis,avg);
}

vfps::meshdata_t vfps::PhaseSpace::variance(const unsigned int axis)
{
	if (axis == 0) {
		projectionToX();
	} else {
		projectionToY();
	}

	if (_moment[axis].size() < 2) {
		_moment[axis].resize(2);
	}

	average(axis);

	meshdata_t avg = _moment[axis][0];
	meshdata_t var = 0;
	for (unsigned int i=0; i<size(axis); i++) {
		var += (meshdata_t(i)-avg)*_projection[axis][i];
	}
	var = var/meshdata_t(size(axis));

	_moment[axis][1] = var;

	return x(axis,var);
}

vfps::meshdata_t* vfps::PhaseSpace::projectionToX() {
	for (unsigned int x=0; x < size(0); x++) {
		_projection[0][x] = 0;
		for (unsigned int y=0; y< size(1); y++) {
			_projection[0][x] += _data[x][y]*_ws[0][y];
		}
	}
	return _projection[0];
}

vfps::meshdata_t* vfps::PhaseSpace::projectionToY() {
	for (unsigned int y=0; y< size(1); y++) {
		_projection[1][y] = 0;
		for (unsigned int x=0; x < size(0); x++) {
			_projection[1][y] += _data[x][y]*_ws[1][x];
		}
	}
	return _projection[1];
}

vfps::PhaseSpace& vfps::PhaseSpace::operator=(const vfps::PhaseSpace& other)
{
	if (_axis[0] == other._axis[0] && _axis[1] == other._axis[1]) {
		std::copy(other._data1D,other._data1D+size(0)*size(1),_data1D);
	} else {
		throw std::range_error("Tried to assign two Mesh2D objects"
							   "with different dimensions.");
	}
	return *this;
}

void vfps::PhaseSpace::renormalize(size_t n, float* args)
{
	float sum=0;
	for( size_t i=0; i<n; i++) {
		sum += args[i];
	}

	if (sum > 0) {
		for( size_t i=0; i<n; i++) {
			args[i] /= sum;
		}
	} else {
		for( size_t i=0; i<n; i++) {
			args[i] = 1.0f/n;
		}
	}
}

void vfps::PhaseSpace::renormalize(size_t n, fixp32* args)
{
	std::vector<fixp32> tmpshare;
	tmpshare.reserve(n);
	fixp32 sum=0;
	for( size_t i=0; i<n; i++) {
		sum += args[i];
		tmpshare.push_back(args[i]);
	}

	if (sum == fixp32(0)) {
		for( size_t i=0; i<n; i++) {
			tmpshare[i] = 1;
		}
		sum = n;
	}
	std::list<std::pair<fixp32,size_t>> remainders;
	fixp32 rensum=0;
	for( size_t i=0; i<n; i++) {
		fixp32 val = tmpshare[i];
		fixp32 remainder = val - sum*(val/sum);
		val /= sum;
		rensum += val;
		args[i] = val;
		remainders.push_back(std::pair<fixp32,size_t>(remainder,i));
	}

	fixp32 rest = fixp32(1) - rensum;
	if (rest > fixp32(0)) {
		// sort remainders descending
		remainders.sort([](	const std::pair<fixp32,size_t> &lhs,
							const std::pair<fixp32,size_t> &rhs)
						-> bool
						{ return lhs.first > rhs.first; }
		);
		do {
			args[remainders.front().second] += std::numeric_limits<fixp32>::epsilon();
			remainders.pop_front();
			rest -= fixp32(std::numeric_limits<fixp32>::epsilon());
		} while (rest > fixp32(0));
	} else  if (rest < fixp32(0)) {
		// sort remainders ascending
		remainders.sort([](	const std::pair<fixp32,size_t> &lhs,
							const std::pair<fixp32,size_t> &rhs)
						-> bool
						{ return lhs.first < rhs.first; }
		);
		do {
			args[remainders.front().second] -= std::numeric_limits<fixp32>::epsilon();
			remainders.pop_front();
			rest += fixp32(std::numeric_limits<fixp32>::epsilon());
		} while (rest > fixp32(0));
	}
}
