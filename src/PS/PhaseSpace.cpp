/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2013-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
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

#include "PS/PhaseSpace.hpp"

#include <numeric>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

vfps::PhaseSpace::PhaseSpace( std::array<meshRuler_ptr, 2> axis
                            , oclhptr_t oclh
                            , const double bunch_charge
                            , const double bunch_current
                            , const uint32_t nbunches
                            , const double zoom
                            , meshdata_t* data
                            )
  : _axis(axis)
  , charge(bunch_charge)
  , current(bunch_current)
  , _nmeshcellsX(nMeshCells(0))
  , _nmeshcellsY(nMeshCells(1))
  , _nbunches(nbunches)
  , _nmeshcells(_nmeshcellsX*_nmeshcellsY)
  , _integralmethod(IntegralMethod::simpson)
  , _bunchpopulation(Array::array1<integral_t>(_nbunches))
  , _integral(0)
  , _projection(Array::array3<projection_t>(2U,_nbunches,_nmeshcellsX))
  , _data(_nbunches,_nmeshcellsX,_nmeshcellsY)
  , _moment(Array::array3<meshaxis_t>(2U,4U,_nbunches))
  , _rms(Array::array2<meshaxis_t>(2U,_nbunches))
  , _ws(simpsonWeights())
  , _oclh(nullptr) // OpenCL will be disabled during dirst initialization steps
  #ifdef INOVESA_USE_OPENGL
  , projectionX_glbuf(0)
  #endif // INOVESA_USE_OPENGL
  #ifdef INOVESA_ENABLE_CLPROFILING
  , xProjEvents(std::make_unique<cl::vector<cl::Event*>>())
  , integEvents(std::make_unique<cl::vector<cl::Event*>>())
  , syncPSEvents(std::make_unique<cl::vector<cl::Event*>>())
  #endif // INOVESA_ENABLE_CLPROFILING
{
    if (data != nullptr) {
        std::copy(data,data+_data.size(),_data());
    } else {
        for (uint32_t n=0; n<_nbunches; n++) {
            gaus(0,n,zoom); // creates gaussian for x axis
            gaus(1,n,zoom); // creates gaussian for y axis
        }

        createFromProjections();
    }

    #ifdef INOVESA_CHG_BUNCH
    std::random_device seed;
    std::default_random_engine engine(seed());

    std::uniform_real_distribution<> xdist(0.0,1.0);
    std::normal_distribution<> ydist(0.0,1.0);


    constexpr meshindex_t nParticles = UINT32_MAX;
    constexpr float amplitude = 2.0f;
    constexpr float pulselen = 1.90e-3f;
    meshindex_t pulsepix = std::ceil(5*pulselen/2.35f/pmax*ps_size);
    constexpr float wavelen = 6.42e-5f;

    meshindex_t x = 0;
    while (x < ps_size/2-pulsepix) {
        for (meshindex_t y = 0; y < ps_size; y++) {
            (*mesh)[x][y]
                =    std::exp(-std::pow((float(x)/ps_size-0.5f)*qmax,2.0f)/2.0f)
                *    std::exp(-std::pow((float(y)/ps_size-0.5f)*pmax,2.0f)/2.0f);
        }
        x++;
    }
    while (x < ps_size/2+pulsepix) {
        meshdata_t weight = std::sqrt(2*pi<meshdata_t>())*ps_size/pmax/nParticles
                * std::exp(-std::pow((float(x)/ps_size-0.5f)*qmax,2.0f)/2.0f);
        for (size_t i=0; i<nParticles; i++) {
            float xf = x+xdist(engine);
            float yf = ydist(engine)
                            + std::exp(-std::pow(xf/(std::sqrt(2)*pulselen/2.35f),2))
                            * amplitude * std::sin(2*pi<meshdata_t>()*xf/wavelen);
            meshindex_t y = std::lround((yf/pmax+0.5f)*ps_size);
            if (y < ps_size) {
                (*mesh)[x][y] += weight;
            }
        }
        x++;
    }
    while (x < ps_size) {
        for (meshindex_t y = 0; y < ps_size; y++) {
            (*mesh)[x][y]
                =    std::exp(-std::pow((float(x)/ps_size-0.5f)*qmax,2.0f)/2.0f)
                *    std::exp(-std::pow((float(y)/ps_size-0.5f)*pmax,2.0f)/2.0f);
        }
        x++;
    }
    #endif // INOVESA_CHG_BUNCH

    _oclh = oclh; // now, OpenCL can be used
    #ifdef INOVESA_USE_OPENCL
    if (_oclh) {
    try {
        data_buf = cl::Buffer(_oclh->context,
                            CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                            sizeof(meshdata_t)*_nbunches*nMeshCells(0)*nMeshCells(1),
                           _data());
        #ifdef INOVESA_USE_OPENGL
        if (_oclh->OpenGLSharing()) {
            glGenBuffers(1, &projectionX_glbuf);
            glBindBuffer(GL_ARRAY_BUFFER,projectionX_glbuf);
            glBufferData( GL_ARRAY_BUFFER
                        , nMeshCells(0)*sizeof(*_projection[0])
                        , 0, GL_DYNAMIC_DRAW);
            projectionX_clbuf = cl::BufferGL( _oclh->context,CL_MEM_READ_WRITE
                                            , projectionX_glbuf);
        } else
        #endif // INOVESA_USE_OPENGL
        {
            projectionX_clbuf = cl::Buffer(_oclh->context,
                                     CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                     sizeof(projection_t)*_nmeshcellsX,
                                     _projection[0]);
        }
        bunchpop_buf = cl::Buffer(_oclh->context,
                                     CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                                     sizeof(integral_t),
                                  &_integral);
        ws_buf = cl::Buffer(_oclh->context,
                            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                            sizeof(meshdata_t)*_nmeshcellsX,
                            _ws());

        _clProgProjX  = _oclh->prepareCLProg(cl_code_projection_x);
        _clKernProjX = cl::Kernel(_clProgProjX, "projectionX");
        _clKernProjX.setArg(0, data_buf);
        _clKernProjX.setArg(1, ws_buf);
        _clKernProjX.setArg(2, _nmeshcellsY);
        _clKernProjX.setArg(3, projectionX_clbuf);

        _clProgIntegral = _oclh->prepareCLProg(cl_code_integral);
        _clKernIntegral = cl::Kernel(_clProgIntegral, "integral");
        _clKernIntegral.setArg(0, projectionX_clbuf);
        _clKernIntegral.setArg(1, ws_buf);
        _clKernIntegral.setArg(2, _nmeshcellsX);
        _clKernIntegral.setArg(3, bunchpop_buf);
    } catch (cl::Error &e) {
        std::cerr << "Error: " << e.what() << std::endl
                  << "Shutting down OpenCL." << std::endl;
        _oclh.reset();
    }
    }
#endif
}

vfps::PhaseSpace::PhaseSpace( meshRuler_ptr axis0
                            , meshRuler_ptr axis1
                            , oclhptr_t oclh
                            , const double bunch_charge
                            , const double bunch_current
                            , uint32_t nbunches
                            , const double zoom
                            , vfps::meshdata_t *data
                            ) :
    PhaseSpace( {{axis0,axis1}}
              , oclh
              , bunch_charge,bunch_current,nbunches,zoom,data)
{
}

vfps::PhaseSpace::PhaseSpace(meshindex_t ps_size
                            , meshaxis_t qmin, meshaxis_t qmax, double qscale
                            , meshaxis_t pmin, meshaxis_t pmax, double pscale
                            , oclhptr_t oclh
                            , const double bunch_charge
                            , const double bunch_current
                            , const uint32_t nbunches
                            , const double zoom, meshdata_t *data
                            )
  : PhaseSpace( meshRuler_ptr(new Ruler<meshaxis_t>(ps_size,qmin,qmax,{{"Meter",qscale}}))
              , meshRuler_ptr(new Ruler<meshaxis_t>(ps_size,pmin,pmax,{{"ElectronVolt",pscale}}))
              , oclh
              , bunch_charge,bunch_current, nbunches, zoom, data)
{}

vfps::PhaseSpace::PhaseSpace(const vfps::PhaseSpace& other) :
    PhaseSpace( other._axis
              , other._oclh
              , other.charge
              , other.current
              , other._nbunches
              , 1
              , other._data
              )
{
}

vfps::PhaseSpace::~PhaseSpace() noexcept
{
    #ifdef INOVESA_ENABLE_CLPROFILING
    if (_oclh) {
        _oclh->saveTimings(xProjEvents.get(),"xProjPS");
        _oclh->saveTimings(integEvents.get(),"integPS");
        _oclh->saveTimings(syncPSEvents.get(),"syncPS");
        for (auto ev : *xProjEvents) {
            delete ev;
        }
        for (auto ev : *integEvents) {
            delete ev;
        }
        for (auto ev : *syncPSEvents) {
            delete ev;
        }
    }
    #endif
}

void vfps::PhaseSpace::integrate()
{
    #ifdef INOVESA_USE_OPENCL
    if (_oclh) {
        _oclh->enqueueNDRangeKernel( _clKernIntegral
                                  , cl::NullRange
                                  , cl::NDRange(1)
                                  #ifdef INOVESA_ENABLE_CLPROFILING
                                  , cl::NullRange
                                  , nullptr
                                  , nullptr
                                  , integEvents.get()
                                  #endif // INOVESA_ENABLE_CLPROFILING
                                  );
    } else
    #endif
    {
    for (uint32_t n=0; n<_nbunches; n++) {
        switch (_integralmethod) {
        case IntegralMethod::sum:
            _bunchpopulation[n] = std::accumulate(_projection[0][n].begin(),
                                           _projection[0][n].end(),
                                           static_cast<integral_t>(0));
            break;
        case IntegralMethod::simpson:
            _bunchpopulation[n] = std::inner_product(_projection[0][n].begin(),
                                              _projection[0][n].end(),
                                              _ws.begin(),
                                              static_cast<integral_t>(0));
            break;
        }
    }
    _integral = std::accumulate( _bunchpopulation.begin()
                               , _bunchpopulation.end()
                               , static_cast<integral_t>(0));
    }
}

void vfps::PhaseSpace::average(const uint_fast8_t axis)
{
    if (axis == 0) {
        #ifdef INOVESA_USE_OPENCL
        if (_oclh) {
        _oclh->enqueueReadBuffer( projectionX_clbuf,CL_TRUE,0
                                , sizeof(projection_t)*nMeshCells(0)
                                , _projection[0]);
        }
        #endif
    }
    for (uint32_t n=0; n<_nbunches; n++) {
        integral_t avg = 0;
        for (size_t i=0; i<nMeshCells(axis); i++) {
            avg += _projection[axis][n][i]*_qp(axis,i);
        }

        // _projection is normalized in p/q coordinates
        avg *= getDelta(axis);

        _moment[axis][0][n] = avg;
    }
}

void vfps::PhaseSpace::variance(const uint_fast8_t axis)
{
    average(axis);
    for (uint32_t n=0; n<_nbunches; n++) {
        meshdata_t var = 0;
        for (size_t i=0; i<nMeshCells(axis); i++) {
            var += _projection[axis][n][i]*std::pow(_qp(axis,i)-_moment[axis][0][n],2);
        }

        // _projection is normalized in p/q coordinates
        var *= getDelta(axis);

        _moment[axis][1][n] = var;
        _rms[axis][n] = std::sqrt(var);
    }
}

void vfps::PhaseSpace::updateXProjection() {
#ifdef INOVESA_USE_OPENCL
    if (_oclh) {
        _oclh->enqueueNDRangeKernel(_clKernProjX
                                  , cl::NullRange
                                  , cl::NDRange(nMeshCells(0))
                                  # ifdef INOVESA_ENABLE_CLPROFILING
                                  , cl::NullRange
                                  , nullptr
                                  , nullptr
                                  , xProjEvents.get()
                                  #endif // INOVESA_ENABLE_CLPROFILING
                                  );
        _oclh->enqueueBarrier();
        #ifdef INOVESA_SYNC_CL
        _oclh->enqueueReadBuffer(projectionX_buf,CL_TRUE,0,
                                      sizeof(projection_t)*nMeshCells(0),
                                      _projection[0].data());
        #endif
    } else
    #endif
    {
    switch (_integralmethod) {
    case IntegralMethod::sum:
        for (size_t n=0; n < _nbunches; n++) {
            for (size_t x=0; x < nMeshCells(0); x++) {
                _projection[0][n][x]
                        = std::accumulate(_data[n][x].begin(),
                                          _data[n][x].end(),
                                          static_cast<integral_t>(0));
                _projection[0][n][x] /= length(1);
            }
        }
        break;
    case IntegralMethod::simpson:
          for (size_t n=0; n < _nbunches; n++) {
              for (size_t x=0; x < nMeshCells(0); x++) {
                  _projection[0][n][x]
                          = std::inner_product(_data[n][x].begin(),
                                               _data[n][x].end(),
                                               _ws.begin(),
                                               static_cast<integral_t>(0));
              }
          }
          break;
    }
    }
}

void vfps::PhaseSpace::updateYProjection() {
    #ifdef INOVESA_USE_OPENCL
    if (_oclh) {
        _oclh->enqueueReadBuffer
            (data_buf,CL_TRUE,0,sizeof(meshdata_t)*nMeshCells(),_data());
    }
    #endif
    switch (_integralmethod) {
    case IntegralMethod::sum:
        for (size_t n=0; n < _nbunches; n++) {
            for (size_t y=0; y< nMeshCells(1); y++) {
                _projection[1][n][y] = 0;

                for (size_t x=0; x< nMeshCells(0); x++) {
                    _projection[1][n][y] += _data[n][x][y];
                }
                _projection[1][n][y] /= length(0);
            }
        }
        break;
    case IntegralMethod::simpson:
        for (size_t n=0; n < _nbunches; n++) {
            for (size_t y=0; y< nMeshCells(1); y++) {
                _projection[1][n][y] = 0;
                for (size_t x=0; x< nMeshCells(0); x++) {
                    _projection[1][n][y] += _data[n][x][y]*_ws[x];
                }
            }
        }
        break;
    }
}

Array::array1<vfps::integral_t> vfps::PhaseSpace::normalize()
{
    integrate();

    #ifdef INOVESA_USE_OPENCL
    syncCLMem(OCLH::clCopyDirection::dev2cpu);
    #endif // INOVESA_USE_OPENCL

    _data /= _integral;

    #ifdef INOVESA_USE_OPENCL
    if (_oclh) {
        _oclh->enqueueWriteBuffer
            (data_buf,CL_TRUE,0,
             sizeof(meshdata_t)*nMeshCells(),_data());
    }
    #endif // INOVESA_USE_OPENCL
    return _bunchpopulation;
}

vfps::PhaseSpace& vfps::PhaseSpace::operator=(vfps::PhaseSpace other)
{
    swap(*this,other);
    return *this;
}

#ifdef INOVESA_USE_OPENCL
void vfps::PhaseSpace::syncCLMem(OCLH::clCopyDirection dir,cl::Event* evt)
{
    if (_oclh) {
    switch (dir) {
    case OCLH::clCopyDirection::cpu2dev:
        _oclh->enqueueWriteBuffer
            (data_buf,CL_TRUE,0,
             sizeof(meshdata_t)*nMeshCells(),_data(),nullptr,evt);
        break;
    case OCLH::clCopyDirection::dev2cpu:
        _oclh->enqueueReadBuffer
            (data_buf,CL_TRUE,0,sizeof(meshdata_t)*nMeshCells(),_data());
        _oclh->enqueueReadBuffer( projectionX_clbuf,CL_TRUE,0
                                , sizeof(projection_t)*nMeshCells(0)
                                , _projection[0],nullptr,evt);
        _oclh->enqueueReadBuffer
            (bunchpop_buf,CL_TRUE,0,sizeof(integral_t),&_bunchpopulation,
            nullptr,evt
            #ifdef INOVESA_ENABLE_CLPROFILING
            , syncPSEvents.get()
            #endif // INOVESA_ENABLE_CLPROFILING
            );
        break;
    }
    }
}
#endif // INOVESA_USE_OPENCL

void vfps::PhaseSpace::createFromProjections()
{
    _data.Activate();
    for (size_t n=0; n < _nbunches; n++) {
        for (meshindex_t x = 0; x < nMeshCells(0); x++) {
            for (meshindex_t y = 0; y < nMeshCells(1); y++) {
                _data[n][x][y] = _projection[0][n][x]*_projection[1][n][y];
            }
        }
    }
    normalize();
}

void vfps::PhaseSpace::gaus(const uint_fast8_t axis,
                            const uint32_t bunch,
                            const double zoom)
{
    const double zoom2=zoom*zoom;
    for(uint32_t i=0; i<nMeshCells(axis); i++){
        _projection[axis][bunch][i]=std::exp((-0.5)*_axis[axis]->at(i)*_axis[axis]->at(i)/zoom2);
    }
}

const Array::array1<vfps::meshdata_t> vfps::PhaseSpace::simpsonWeights()
{
    Array::array1<vfps::meshdata_t> rv(nMeshCells(0));
    const integral_t ca = 3.;
    integral_t dc = 1;

    const integral_t h03 = getDelta(0)/integral_t(3);
    rv[0] = h03;
    for (size_t x=1; x< nMeshCells(0)-1; x++){
        _ws[x] = h03 * (ca+dc);
        dc = -dc;
    }
    rv[nMeshCells(0)-1] = h03;

    return rv;
}

void vfps::swap(vfps::PhaseSpace& first, vfps::PhaseSpace& second) noexcept
{
    std::swap(first._data, second._data);
}

#ifdef INOVESA_USE_OPENCL
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
#endif // INOVESA_USE_OPENCL

