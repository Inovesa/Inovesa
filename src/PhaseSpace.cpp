/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   *
 * Copyright (c) 2014-2016: Patrik Schönfeldt                                 *
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

vfps::PhaseSpace::PhaseSpace(std::array<Ruler<meshaxis_t>,2> axis,
                             const double Fk,meshindex_t xoffset) :
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

    if (Fk >= 0) {
        if (Fk > 0) {
            haissinski(0,Fk); // 50 iterations haissinski for y axis
        } else {
            haissinski(0,0); // creates gaussian for x axis
        }
        haissinski(1,0); // creates gaussian for y axis

        for (meshindex_t x = 0; x < nMeshCells(0); x++) {
            for (meshindex_t y = 0; y < nMeshCells(1); y++) {
                _data[x][y] = _projection[0][(x+xoffset)%nMeshCells(0)]*_projection[1][y];
            }
        }
    }

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

        _clProgProjX  = OCLH::prepareCLProg(cl_code_projection_x);
        _clKernProjX = cl::Kernel(_clProgProjX, "projectionX");
        _clKernProjX.setArg(0, data_buf);
        _clKernProjX.setArg(1, ws_buf);
        _clKernProjX.setArg(2, _nmeshcellsY);
        _clKernProjX.setArg(3, projectionX_buf);

        _clProgIntegral = OCLH::prepareCLProg(cl_code_integral);
        _clKernIntegral = cl::Kernel(_clProgIntegral, "integral");
        _clKernIntegral.setArg(0, projectionX_buf);
        _clKernIntegral.setArg(1, ws_buf);
        _clKernIntegral.setArg(2, _nmeshcellsX);
        _clKernIntegral.setArg(3, integral_buf);
    }
    #endif
}


vfps::PhaseSpace::PhaseSpace(Ruler<meshaxis_t> axis1, Ruler<meshaxis_t> axis2,
                             const double Fk) :
    PhaseSpace(std::array<Ruler<meshaxis_t>,2>{{axis1,axis2}},Fk)
{}

vfps::PhaseSpace::PhaseSpace(meshindex_t ps_size,
                             meshaxis_t xmin, meshaxis_t xmax,
                             meshaxis_t ymin, meshaxis_t ymax,
                             double xscale, double yscale,
                             const double Fk) :
    PhaseSpace(Ruler<meshaxis_t>(ps_size,xmin,xmax,xscale),
               Ruler<meshaxis_t>(ps_size,ymin,ymax,yscale),
               Fk)
{}

vfps::PhaseSpace::PhaseSpace(const vfps::PhaseSpace& other) :
    PhaseSpace(other._axis,-1)
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
        OCLH::queue.enqueueReadBuffer
            (integral_buf,CL_TRUE,0,sizeof(integral_t),&_integral);
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

vfps::meshaxis_t vfps::PhaseSpace::average(const uint_fast8_t axis)
{
    if (axis == 0) {
        updateXProjection();
    } else {
        updateYProjection();
    }
    integral_t avg = 0;
    for (size_t i=0; i<nMeshCells(axis); i++) {
        avg += _projection[axis][i]*x(axis,i);
    }

    // _projection is normalized in p/q coordinates
    avg *= getDelta(axis);

    _moment[axis][0] = avg;

    return static_cast<meshaxis_t>(avg);
}

vfps::meshdata_t vfps::PhaseSpace::variance(const uint_fast8_t axis)
{
    meshdata_t avg = average(axis);;
    meshdata_t var = 0;
    for (size_t i=0; i<nMeshCells(axis); i++) {
        var += _projection[axis][i]*std::pow(x(axis,i)-avg,2);
    }

    // _projection is normalized in p/q coordinates
    var *= getDelta(axis);

    _moment[axis][1] = var;

    return var;
}

void vfps::PhaseSpace::updateXProjection() {
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        OCLH::queue.enqueueNDRangeKernel (
                    _clKernProjX,
                    cl::NullRange,
                    cl::NDRange(nMeshCells(0)));
        OCLH::queue.enqueueBarrierWithWaitList();
        #ifdef INOVESA_SYNC_CL
        syncCLMem(clCopyDirection::dev2cpu);
        #endif
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

vfps::integral_t vfps::PhaseSpace::normalize()
{
    integral();
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        OCLH::queue.enqueueReadBuffer
            (data_buf,CL_TRUE,0,sizeof(meshdata_t)*nMeshCells(),_data1D);
    }
    #endif // INOVESA_USE_CL
    for (meshindex_t i = 0; i < _nmeshcells; i++) {
        _data1D[i] /= _integral;
    }
    #ifdef INOVESA_USE_CL
    if (OCLH::active) {
        OCLH::queue.enqueueWriteBuffer
            (data_buf,CL_TRUE,0,
             sizeof(meshdata_t)*nMeshCells(),_data1D);
    }
    #endif // INOVESA_USE_CL
    return _integral;
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

void vfps::PhaseSpace::haissinski(const uint_fast8_t x,
                                  const projection_t Fk)
{
    constexpr uint32_t maxloops = 200;
    constexpr double kappamax=0.291030514208;
    constexpr double kappamin=0.01;
    constexpr double p=1.50088;
    constexpr double q=1.05341;
    projection_t kappa = std::max(kappamax*(1.0-std::exp(-p*std::pow(Fk,q))),
                                  kappamin);
    projection_t* I = new projection_t[nMeshCells(x)];
    std::fill_n(I,nMeshCells(x),0);
    projection_t F=-1;
    for(uint32_t k=0;k<=maxloops && F<Fk;k++){ // fuehrt Iterationen durch
        F=0;
        for(uint32_t i=0; i<nMeshCells(x); i++){
            data_t tv=0; // vorheriges t
            _projection[x][i]=kappa*std::exp((-0.5)*_axis[x][i]*_axis[x][i]+I[i]);
            F+=_projection[x][i]*getDelta(x);
            I[i]=0;
            for(uint32_t j=0; j<=i; j++){ //Berechnet neuen Wert I[i]
                data_t tn =std::pow(((j+1)*getDelta(x)),2./3.); // naechstes t
                data_t b  = (tn-tv)*0.5;
                I[i]+= _projection[x][i-j]*b;
                tv =std::pow((j*getDelta(x)),2./3.);
            }
            I[i]*=1.5;
        }
    }

    for (uint32_t i=0;i<nMeshCells(x);i++){ // Normalize Haissinski distribution
        _projection[x][i]/=F;
    }
    delete [] I;
}
#endif // INOVESA_USE_CL

void vfps::swap(vfps::PhaseSpace& first, vfps::PhaseSpace& second) noexcept
{
    std::swap(first._data, second._data);
    std::swap(first._data1D,second._data1D);
}

std::string vfps::PhaseSpace::cl_code_integral = R"(
    __kernel void integral(const __global float* proj,
                           const __global float* ws,
                           const uint xsize,
                           __global float* result)
    {
        float value = 0;
        for (uint x=0; x< xsize; x++) {
            value += proj[x]*ws[x];
        }
        *result = value;
    }
    )";

std::string vfps::PhaseSpace::cl_code_projection_x = R"(
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

