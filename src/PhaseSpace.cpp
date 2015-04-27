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

#include "PhaseSpace.hpp"

vfps::PhaseSpace::PhaseSpace(std::array<Ruler<meshaxis_t>,2> axis) :
	_axis(axis),
	_data1D(new meshdata_t[nMeshCells(0)*nMeshCells(1)]())
{
	_data = new meshdata_t*[nMeshCells(0)];
	for (unsigned int i=0; i<nMeshCells(0); i++) {
		_data[i] = &(_data1D[i*nMeshCells(1)]);
	}
	_projection[0] = new integral_t[nMeshCells(0)];
	_projection[1] = new integral_t[nMeshCells(1)];
	_ws[0] = new meshdata_t[nMeshCells(0)];
	_ws[1] = new meshdata_t[nMeshCells(1)];

	const integral_t ca = 3.;
	integral_t dc = 1;

	const integral_t h03 = getDelta(0)/3.;
	_ws[0][0] = h03;
	for (unsigned int x=1; x< nMeshCells(0)-1; x++){
		_ws[0][x] = h03 * (ca+dc);
		dc = -dc;
	}
	_ws[0][nMeshCells(0)-1] = h03;


	const integral_t h13 = getDelta(1)/3.;
	_ws[1][0] = h13;
	for (unsigned int x=1; x< nMeshCells(1)-1; x++){
		_ws[1][x] = h13 * (ca+dc);
		dc = -dc;
	}
	_ws[1][nMeshCells(0)-1] = h13;
	#ifdef INOVESA_USE_CL
	data_buf = cl::Buffer(OCLH::context,
							CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
							sizeof(meshdata_t)*nMeshCells(0)*nMeshCells(1),
						   _data1D);
	#endif
}


vfps::PhaseSpace::PhaseSpace(Ruler<meshaxis_t> axis1, Ruler<meshaxis_t> axis2) :
	PhaseSpace(std::array<Ruler<meshaxis_t>,2>{{axis1,axis2}})
{}

vfps::PhaseSpace::PhaseSpace(unsigned int ps_size,
							 meshaxis_t xmin, meshaxis_t xmax,
							 meshaxis_t ymin, meshaxis_t ymax) :
	PhaseSpace(Ruler<meshaxis_t>(ps_size,xmin,xmax),
			   Ruler<meshaxis_t>(ps_size,ymin,ymax))
{}

vfps::PhaseSpace::PhaseSpace(const vfps::PhaseSpace& other) :
	PhaseSpace(other._axis)
{ std::copy_n(other._data1D,nMeshCells(0)*nMeshCells(1),_data1D); }

vfps::PhaseSpace::~PhaseSpace()
{
	delete [] _data;
	delete [] _data1D;
	delete [] _projection[0];
	delete [] _projection[1];
	delete [] _ws[0];
	delete [] _ws[1];
}

/*
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
*/

vfps::integral_t* vfps::PhaseSpace::projectionToX() {
	for (unsigned int x=0; x < nMeshCells(0); x++) {
		_projection[0][x] = 0;

		for (unsigned int y=0; y< nMeshCells(1); y++) {
			#if INTEGRAL_TYPE == 1
				_projection[0][x] += _data[x][y];
			#elif INTEGRAL_TYPE == 2
				_projection[0][x] += _data[x][y]*_ws[0][y];
			#endif
		}
		#if INTEGRAL_TYPE == 1
			_projection[0][x] /= size(0);
		#endif
	}
	return _projection[0];
}

vfps::integral_t* vfps::PhaseSpace::projectionToY() {
	for (unsigned int y=0; y< nMeshCells(1); y++) {
		_projection[1][y] = 0;

		for (unsigned int x=0; x< nMeshCells(0); y++) {
			#if INTEGRAL_TYPE == 1
				_projection[1][y] += _data[x][y];
			#elif INTEGRAL_TYPE == 2
				_projection[1][y] += _data[x][y]*_ws[1][x];
			#endif
		}
		#if INTEGRAL_TYPE == 1
			_projection[1][y] /= size(1);
		#endif
	}
	return _projection[1];
}

vfps::PhaseSpace& vfps::PhaseSpace::operator=(vfps::PhaseSpace other)
{
	swap(*this,other);
	return *this;
}

#ifdef INOVESA_USE_CL
void vfps::PhaseSpace::syncCLMem(clCopyDirection dir)
{
	switch (dir) {
	case clCopyDirection::cpu2dev:
		OCLH::queue.enqueueWriteBuffer
			(data_buf,CL_TRUE,0,
			 sizeof(meshdata_t)*nMeshCells(),_data1D);
		break;
	case clCopyDirection::dev2cpu:
		OCLH::queue.enqueueReadBuffer
			(data_buf,CL_TRUE,0,sizeof(meshdata_t)*nMeshCells(),_data1D);
		break;
	}
}
#endif // INOVESA_USE_CL

void vfps::swap(vfps::PhaseSpace& first, vfps::PhaseSpace& second) noexcept
{
	std::swap(first._data, second._data);
	std::swap(first._data1D,second._data1D);
}
