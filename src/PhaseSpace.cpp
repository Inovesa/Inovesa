/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   *
 * Copyright (c) 2014-2015: Patrik Schönfeldt                                 *
 *                                                                            *
 * This file is part of Inovesa.                                              *
 * Inovesa is free software: you can redistribute it and/or modify            *
 * it under the terms of the GNU General Public License as published by       *
 * the Free Software Foundation, either version 3 of the License, or          *
 * (at your option) any later version.                                        *
 *                                                                            *
 * Inovesa is distributed in the hope that it will be useful,                 *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of             *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
 * GNU General Public License for more details.                               *
 *                                                                            *
 * You should have received a copy of the GNU General Public License          *
 * along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           *
 ******************************************************************************/

#include "PhaseSpace.hpp"

vfps::PhaseSpace::PhaseSpace(std::array<Ruler<meshaxis_t>,2> axis) :
    _axis(axis),
    _integral(0),
    _projection(std::array<integral_t*,2> {{ new integral_t[nMeshCells(0)],
                                             new integral_t[nMeshCells(1)]
                                          }}),
    _nmeshcellsX(nMeshCells(0)),
    _nmeshcellsY(nMeshCells(1)),
    _nmeshcells(_nmeshcellsX*_nmeshcellsY),
    _integraltype(IntegralType::simpson),
    _data1D(new meshdata_t[_nmeshcells]())
{
    _data = new meshdata_t*[nMeshCells(0)];
    for (size_t i=0; i<nMeshCells(0); i++) {
        _data[i] = &(_data1D[i*nMeshCells(1)]);
    }
    _ws = new meshdata_t[nMeshCells(0)];

    const integral_t ca = 3.;
    integral_t dc = 1;

    const integral_t h03 = getDelta(0)/integral_t(3);
    _ws[0] = h03;
    for (size_t x=1; x< nMeshCells(0)-1; x++){
        _ws[x] = h03 * (ca+dc);
        dc = -dc;
    }
    _ws[nMeshCells(0)-1] = h03;

    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        data_buf = cl::Buffer(OCLH::context,
                            CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                            sizeof(meshdata_t)*nMeshCells(0)*nMeshCells(1),
                           _data1D);
        projectionX_buf = cl::Buffer(OCLH::context,
                                     CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                     sizeof(projection_t)*_nmeshcellsX,
                                     _projection[0]);
        integral_buf = cl::Buffer(OCLH::context,
                                     CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                     sizeof(integral_t),&_integral);
        ws_buf = cl::Buffer(OCLH::context,
                            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                            sizeof(meshdata_t)*_nmeshcellsX,
                            _ws);
        std::string cl_code_projection_x = R"(
            __kernel void projectionX(const __global float* mesh,
                                      const __global float* ws,
                                      const uint ysize,
                                      __global float* proj)
            {
                float value = 0;
                const uint x = get_global_id(0);
                const uint meshoffs = x*ysize;
                for (uint y=0; y< ysize; y++) {
                    value += mesh[meshoffs+y]*ws[y];
                }
                proj[x] = value;
            }
            )";
        _clProgProjX  = OCLH::prepareCLProg(cl_code_projection_x);
        _clKernProjX = cl::Kernel(_clProgProjX, "projectionX");
        _clKernProjX.setArg(0, data_buf);
        _clKernProjX.setArg(1, ws_buf);
        _clKernProjX.setArg(2, _nmeshcellsY);
        _clKernProjX.setArg(3, projectionX_buf);

        std::string cl_code_integral = R"(
            __kernel void integral(const __global float* proj,
                                   const __global float* ws,
                                   const uint xsize,
                                   float result)
            {
                float value = 0;
                for (uint x=0; x< xsize; x++) {
                    value += proj[x]*ws[x];
                }
                result = value;
            }
            )";
        _clProgIntegral = OCLH::prepareCLProg(cl_code_integral);
        _clKernIntegral = cl::Kernel(_clProgIntegral, "integral");
        _clKernIntegral.setArg(0, projectionX_buf);
        _clKernIntegral.setArg(1, ws_buf);
        _clKernIntegral.setArg(2, _nmeshcellsX);
        _clKernIntegral.setArg(3, _integral);
    }
    #endif
}


vfps::PhaseSpace::PhaseSpace(Ruler<meshaxis_t> axis1, Ruler<meshaxis_t> axis2) :
    PhaseSpace(std::array<Ruler<meshaxis_t>,2>{{axis1,axis2}})
{}

vfps::PhaseSpace::PhaseSpace(meshindex_t ps_size,
                             meshaxis_t xmin, meshaxis_t xmax,
                             meshaxis_t ymin, meshaxis_t ymax,
                             double xscale, double yscale) :
    PhaseSpace(Ruler<meshaxis_t>(ps_size,xmin,xmax,xscale),
               Ruler<meshaxis_t>(ps_size,ymin,ymax,yscale))
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
    delete [] _ws;
}

vfps::integral_t vfps::PhaseSpace::integral()
{
    updateXProjection();
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        OCLH::queue.enqueueNDRangeKernel (
                    _clKernIntegral,
                    cl::NullRange,
                    cl::NDRange(1));
    } else
    #endif
    {
    _integral = 0;
    switch (_integraltype) {
    case IntegralType::sum:
        for (size_t x=0; x< nMeshCells(0); x++) {
            _integral += _projection[0][x];
        }
        break;
    case IntegralType::simpson:
        for (size_t x=0; x< nMeshCells(0); x++) {
            _integral += _projection[0][x]*static_cast<integral_t>(_ws[x]);
        }
        break;
    }
    }
    return _integral;
}

/*
vfps::meshdata_t vfps::PhaseSpace::average(const uint_fast8_t axis)
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
    for (size_t i=0; i<size(axis); i++) {
        avg += meshdata_t(i)*_projection[axis][i];
    }
    avg = avg/meshdata_t(size(axis));

    _moment[axis][0] = avg;

    return x(axis,avg);
}

vfps::meshdata_t vfps::PhaseSpace::variance(const uint_fast8_t axis)
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
    for (size_t i=0; i<size(axis); i++) {
        var += (meshdata_t(i)-avg)*_projection[axis][i];
    }
    var = var/meshdata_t(size(axis));

    _moment[axis][1] = var;

    return x(axis,var);
}
*/

void vfps::PhaseSpace::updateXProjection(bool sync) {
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        OCLH::queue.enqueueNDRangeKernel (
                    _clKernProjX,
                    cl::NullRange,
                    cl::NDRange(nMeshCells(0)));
        OCLH::queue.enqueueBarrierWithWaitList();
        if (sync) {
            OCLH::queue.enqueueReadBuffer (projectionX_buf,CL_TRUE,0,
                                           sizeof(projection_t)*nMeshCells(0),
                                           _projection[0]);
        }
    } else
    #endif
    {
    switch (_integraltype) {
    case IntegralType::sum:
        for (size_t x=0; x < nMeshCells(0); x++) {
            _projection[0][x] = 0;
            for (size_t y=0; y< nMeshCells(1); y++) {
                _projection[0][x] += _data[x][y];
            }
            _projection[0][x] /= size(1);
        }
        break;
    case IntegralType::simpson:
          for (size_t x=0; x < nMeshCells(0); x++) {
              _projection[0][x] = 0;
              for (size_t y=0; y< nMeshCells(1); y++) {
                _projection[0][x] += _data[x][y]*_ws[y];
              }
          }
          break;
        }
    }
}

void vfps::PhaseSpace::updateYProjection() {
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        syncCLMem(clCopyDirection::dev2cpu);
    }
    #endif
    switch (_integraltype) {
    case IntegralType::sum:
        for (size_t y=0; y< nMeshCells(1); y++) {
            _projection[1][y] = 0;

            for (size_t x=0; x< nMeshCells(0); x++) {
                _projection[1][y] += _data[x][y];
            }
            _projection[1][y] /= size(0);
        }
        break;
    case IntegralType::simpson:
        for (size_t y=0; y< nMeshCells(1); y++) {
            _projection[1][y] = 0;
            for (size_t x=0; x< nMeshCells(0); x++) {
                _projection[1][y] += _data[x][y]*_ws[x];
            }
        }
        break;
    }
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
        OCLH::queue.enqueueReadBuffer(projectionX_buf,CL_TRUE,0,
                                      sizeof(projection_t)*nMeshCells(0),
                                      _projection[0]);
        break;
    }
}
#endif // INOVESA_USE_CL

void vfps::swap(vfps::PhaseSpace& first, vfps::PhaseSpace& second) noexcept
{
    std::swap(first._data, second._data);
    std::swap(first._data1D,second._data1D);
}
