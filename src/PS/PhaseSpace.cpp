// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Tobias Boltz
 * Copyright (c) Karlsruhe Institute of Technology
 */

#include "PS/PhaseSpace.hpp"

#include <iostream>
#include <numeric>
#include <stdexcept>

#include <boost/math/constants/constants.hpp>
using boost::math::constants::pi;

vfps::PhaseSpace::PhaseSpace( std::array<meshRuler_ptr, 2> axis
                            , oclhptr_t oclh
                            , const double beam_charge
                            , const double beam_current
                            , const std::vector<integral_t> filling
                            , const double zoom
                            , const meshdata_t* data
                            )
  : _axis(axis)
  , charge(beam_charge)
  , current(beam_current)
  , _filling_set(filling.begin(),filling.end())
  , _filling(std::vector<integral_t>(_nbunches))
  , _integral(1)
  , _projection(boost::extents[2U][_nbunches][_nmeshcellsX])
  , _data(boost::extents[_nbunches][_nmeshcellsX][_nmeshcellsY])
  , _moment(boost::extents[2U][4U][_nbunches])
  , _rms(boost::extents[2U][_nbunches])
  , _ws(simpsonWeights())
  , _oclh(nullptr) // OpenCL will be disabled during dirst initialization steps
  #if INOVESA_USE_OPENGL == 1
  , projectionX_glbuf(0)
  #endif // INOVESA_USE_OPENGL
  #if INOVESA_ENABLE_CLPROFILING == 1
  , xProjEvents(std::make_unique<cl::vector<cl::Event*>>())
  , integEvents(std::make_unique<cl::vector<cl::Event*>>())
  , syncPSEvents(std::make_unique<cl::vector<cl::Event*>>())
  #endif // INOVESA_ENABLE_CLPROFILING
{

    if (std::round(1e5*std::accumulate( filling.begin(), filling.end(),
                                        static_cast<integral_t>(0)))
            != 1e5) {
        throw std::invalid_argument("Argument \"filling\" not normalized.");
    }
    if (data != nullptr) {
        std::copy(data,data+_totalmeshcells,_data.data());
    } else {
        for (meshindex_t  n=0; n<_nbunches; n++) {
            gaus(0,n,zoom); // creates gaussian for x axis
            gaus(1,n,zoom); // creates gaussian for y axis
        }

        createFromProjections();
    }

    updateXProjection();
    updateYProjection();
    integrate();

    _oclh = oclh; // now, OpenCL can be used
    #if INOVESA_USE_OPENCL == 1
    if (_oclh) {
    try {
        data_buf = cl::Buffer(_oclh->context,
                            CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                            sizeof(meshdata_t)*_totalmeshcells,
                           _data.data());
        #if INOVESA_USE_OPENGL == 1
        if (_oclh->OpenGLSharing()) {
            glGenBuffers(1, &projectionX_glbuf);
            glBindBuffer(GL_ARRAY_BUFFER,projectionX_glbuf);
            glBufferData( GL_ARRAY_BUFFER
                        , _nmeshcellsX*sizeof(decltype(_projection)::value_type)
                        , 0, GL_DYNAMIC_DRAW);
            projectionX_clbuf = cl::BufferGL( _oclh->context,CL_MEM_READ_WRITE
                                            , projectionX_glbuf);
        } else
        #endif // INOVESA_USE_OPENGL
        {
            // _projection is used for _projection[0]
            projectionX_clbuf = cl::Buffer(
                        _oclh->context,
                        CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                        _nmeshcellsX*sizeof(decltype(_projection)::value_type),
                        _projection.data());
        }
        filling_buf = cl::Buffer( _oclh->context
                                 , CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR
                                 , _nbunches*sizeof(integral_t)
                                 , _filling.data());
        ws_buf = cl::Buffer(_oclh->context,
                            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                            sizeof(decltype(_ws)::value_type)*_nmeshcellsX,
                            const_cast<decltype(_ws)::value_type*>(_ws.data()));

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
        _clKernIntegral.setArg(3, filling_buf);
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
                            , const double beam_charge
                            , const double beam_current
                            , const std::vector<integral_t> filling
                            , const double zoom
                            , const meshdata_t* data
                            ) :
    PhaseSpace( {{axis0,axis1}}
              , oclh
              , beam_charge,beam_current,filling,zoom,data)
{
}

vfps::PhaseSpace::PhaseSpace( meshaxis_t qmin, meshaxis_t qmax, double qscale
                            , meshaxis_t pmin, meshaxis_t pmax, double pscale
                            , oclhptr_t oclh
                            , const double beam_charge
                            , const double beam_current
                            , const std::vector<integral_t> filling
                            , const double zoom, const meshdata_t* data
                            )
  : PhaseSpace( meshRuler_ptr(new Ruler<meshaxis_t>(PhaseSpace::nx,qmin,qmax,
                {{"Meter",qscale}}))
              , meshRuler_ptr(new Ruler<meshaxis_t>(PhaseSpace::ny,pmin,pmax,
                {{"ElectronVolt",pscale}}))
              , oclh
              , beam_charge,beam_current, filling, zoom, data)
{}

vfps::PhaseSpace::PhaseSpace(const vfps::PhaseSpace& other) :
    PhaseSpace( other._axis
              , other._oclh
              , other.charge
              , other.current
              , other._filling_set
              , 1 // zoom
              , other._data.data()
              )
{
}

vfps::PhaseSpace::~PhaseSpace() noexcept
{
    #if INOVESA_ENABLE_CLPROFILING == 1
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
    #if INOVESA_USE_OPENCL == 1
    if (_oclh) {
        _oclh->enqueueNDRangeKernel( _clKernIntegral
                                  , cl::NullRange
                                  , cl::NDRange(1)
                                  #if INOVESA_ENABLE_CLPROFILING == 1
                                  , cl::NullRange
                                  , nullptr
                                  , nullptr
                                  , integEvents.get()
                                  #endif // INOVESA_ENABLE_CLPROFILING
                                  );
    } else
    #endif
    {
    for (meshindex_t n=0; n<_nbunches; n++) {
        _filling[n] = std::inner_product(_projection[0][n].begin(),
                                         _projection[0][n].end(),
                                         _ws.begin(),
                                         static_cast<integral_t>(0));
    }
    _integral = std::accumulate( _filling.begin()
                               , _filling.end()
                               , static_cast<integral_t>(0));
    }
}

void vfps::PhaseSpace::average(const uint_fast8_t axis)
{
    if (axis == 0) {
        #if INOVESA_USE_OPENCL == 1
        if (_oclh) {
        // _projection is used for _projection[0]
        _oclh->enqueueReadBuffer( projectionX_clbuf,CL_TRUE,0
                                , sizeof(projection_t)*_nmeshcellsX
                                , _projection.data());
        }
        #endif
    }

    const meshindex_t maxi = (axis==0)? _nmeshcellsX : _nmeshcellsY;
    for (meshindex_t n=0; n<_nbunches; n++) {
        integral_t avg = 0;
        if (_filling_set[n] > 0) {
            for (size_t i=0; i<maxi; i++) {
                avg += _projection[axis][n][i]*_qp(axis,i);
            }

            // _projection is normalized in p/q coordinates with respect to charge
            avg *= getDelta(axis)/_filling[n];
        }

        _moment[axis][0][n] = avg;
    }
}

void vfps::PhaseSpace::variance(const uint_fast8_t axis)
{
    average(axis);
    const meshindex_t maxi = (axis==0)? _nmeshcellsX : _nmeshcellsY;
    for (meshindex_t n=0; n<_nbunches; n++) {
        meshdata_t var = 0;
        if (_filling_set[n] > 0) {
            for (size_t i=0; i<maxi; i++) {
                var += _projection[axis][n][i]*std::pow(_qp(axis,i)
                                                        -_moment[axis][0][n],2);
            }

            // _projection is normalized in p/q coordinates with respect to charge
            var *= getDelta(axis)/_filling[n];
        }

        _moment[axis][1][n] = var;
        _rms[axis][n] = std::sqrt(var);
    }
}

void vfps::PhaseSpace::updateXProjection() {
#if INOVESA_USE_OPENCL == 1
    if (_oclh) {
        _oclh->enqueueNDRangeKernel(_clKernProjX
                                  , cl::NullRange
                                  , cl::NDRange(_nmeshcellsX)
                                  # if INOVESA_ENABLE_CLPROFILING == 1
                                  , cl::NullRange
                                  , nullptr
                                  , nullptr
                                  , xProjEvents.get()
                                  #endif // INOVESA_ENABLE_CLPROFILING
                                  );
        _oclh->enqueueBarrier();
        #if INOVESA_SYNC_CL == 1
        _oclh->enqueueReadBuffer(projectionX_buf,CL_TRUE,0,
                                      sizeof(projection_t)*_nmeshcellsX,
                                      _projection[0].data());
        #endif
    } else
    #endif
    {
        for (size_t n=0; n < _nbunches; n++) {
            for (size_t x=0; x < _nmeshcellsX; x++) {
                _projection[0][n][x]
                      = std::inner_product(_data[n][x].begin(),
                                           _data[n][x].end(),
                                           _ws.begin(),
                                           static_cast<integral_t>(0));
            }
        }
    }
}

void vfps::PhaseSpace::updateYProjection() {
    #if INOVESA_USE_OPENCL == 1
    if (_oclh) {
        _oclh->enqueueReadBuffer
            (data_buf,CL_TRUE,0,sizeof(meshdata_t)*_nmeshcells,_data.data());
    }
    #endif
    for (size_t n=0; n < _nbunches; n++) {
        for (size_t y=0; y< _nmeshcellsY; y++) {
            _projection[1][n][y] = 0;
            for (size_t x=0; x< _nmeshcellsX; x++) {
                _projection[1][n][y] += _data[n][x][y]*_ws[x];
            }
        }
    }
}

const std::vector<vfps::integral_t>& vfps::PhaseSpace::normalize()
{
    #if INOVESA_USE_OPENCL == 1
    syncCLMem(OCLH::clCopyDirection::dev2cpu);
    #endif // INOVESA_USE_OPENCL

    for (size_t n=0; n < _nbunches; n++) {
        if (_filling_set[n] > 0) {
            for (meshindex_t x = 0; x < _nmeshcellsX; x++) {
                for (meshindex_t y = 0; y < _nmeshcellsY; y++) {
                    _data[n][x][y] *= _filling_set[n]/_filling[n];
                }
            }
        } else {
            for (meshindex_t x = 0; x < _nmeshcellsX; x++) {
                for (meshindex_t y = 0; y < _nmeshcellsY; y++) {
                    _data[n][x][y] = 0;
                }
            }
        }
    }

    #if INOVESA_USE_OPENCL == 1
    if (_oclh) {
        _oclh->enqueueWriteBuffer
            (data_buf,CL_TRUE,0,
             sizeof(meshdata_t)*_nmeshcells,_data.data());
    }
    #endif // INOVESA_USE_OPENCL
    return _filling;
}

vfps::PhaseSpace& vfps::PhaseSpace::operator =(vfps::PhaseSpace other)
{
    other.swap(*this);
    return *this;
}

#if INOVESA_USE_OPENCL == 1
void vfps::PhaseSpace::syncCLMem(OCLH::clCopyDirection dir,cl::Event* evt)
{
    if (_oclh) {
    switch (dir) {
    case OCLH::clCopyDirection::cpu2dev:
        _oclh->enqueueWriteBuffer
            (data_buf,CL_TRUE,0,
             sizeof(meshdata_t)*_nmeshcells,_data.data(),nullptr,evt);
        break;
    case OCLH::clCopyDirection::dev2cpu:
        _oclh->enqueueReadBuffer
            (data_buf,CL_TRUE,0,sizeof(meshdata_t)*_nmeshcells,_data.data());
        // _projection is used for _projection[0]
        _oclh->enqueueReadBuffer( projectionX_clbuf,CL_TRUE,0
                                , sizeof(projection_t)*_nmeshcellsX
                                , _projection.data(),nullptr,evt);
        _oclh->enqueueReadBuffer
            (filling_buf,CL_TRUE,0,sizeof(integral_t),&_filling,
            nullptr,evt
            #if INOVESA_ENABLE_CLPROFILING == 1
            , syncPSEvents.get()
            #endif // INOVESA_ENABLE_CLPROFILING
            );
        break;
    }
    }
}
#endif // INOVESA_USE_OPENCL

const vfps::meshindex_t& vfps::PhaseSpace::nx = vfps::PhaseSpace::_nmeshcellsX;

const vfps::meshindex_t& vfps::PhaseSpace::ny = vfps::PhaseSpace::_nmeshcellsY;

const vfps::meshindex_t& vfps::PhaseSpace::nb = vfps::PhaseSpace::_nbunches;

const vfps::meshindex_t& vfps::PhaseSpace::nxy = vfps::PhaseSpace::_nmeshcells;

const vfps::meshindex_t& vfps::PhaseSpace::nxyb
    = vfps::PhaseSpace::_totalmeshcells;

void vfps::PhaseSpace::setSize(const meshindex_t x,
                               const meshindex_t b )
{
    if (vfps::PhaseSpace::_firstinit) {
        vfps::PhaseSpace::_firstinit = false;
        vfps::PhaseSpace::_nmeshcellsX = x;
        vfps::PhaseSpace::_nmeshcellsY = x;
        vfps::PhaseSpace::_nbunches = b;
        vfps::PhaseSpace::_nmeshcells = x*x;
        vfps::PhaseSpace::_totalmeshcells = x*x*b;
    } else {
        std::cerr << "vfps::PhaseSpace::setSize() "
                     "is meant to be called only once." << std::endl;
    }
}

bool vfps::PhaseSpace::_firstinit(true);

vfps::meshindex_t vfps::PhaseSpace::_nmeshcellsX(0);

vfps::meshindex_t vfps::PhaseSpace::_nmeshcellsY(0);

vfps::meshindex_t vfps::PhaseSpace::_nbunches(0);

vfps::meshindex_t vfps::PhaseSpace::_nmeshcells(0);

vfps::meshindex_t vfps::PhaseSpace::_totalmeshcells(0);

void vfps::PhaseSpace::createFromProjections()
{
    for (size_t n=0; n < _nbunches; n++) {
        for (meshindex_t x = 0; x < _nmeshcellsX; x++) {
            for (meshindex_t y = 0; y < _nmeshcellsY; y++) {
                _data[n][x][y] = _projection[0][n][x]*_projection[1][n][y];
            }
        }
    }
    integrate();
    normalize();
}

void vfps::PhaseSpace::gaus(const uint_fast8_t axis,
                            const meshindex_t bunch,
                            const double zoom)
{
    const double zoom2=zoom*zoom;
    const meshindex_t maxi = (axis==0)? _nmeshcellsX : _nmeshcellsY;
    for(meshindex_t i=0; i<maxi; i++){
        _projection[axis][bunch][i]
                = boost::math::constants::one_div_root_two_pi<double>()
                * std::exp((-0.5)*_axis[axis]->at(i)*_axis[axis]->at(i)/zoom2);
    }
}

const std::vector<vfps::meshdata_t> vfps::PhaseSpace::simpsonWeights()
{
    std::vector<vfps::meshdata_t> rv(_nmeshcellsX);
    const integral_t ca = 3.;
    integral_t dc = 1;

    const integral_t h03 = getDelta(0)/integral_t(3);
    rv[0] = h03;
    for (size_t x=1; x< _nmeshcellsX-1; x++){
        rv[x] = h03 * (ca+dc);
        dc = -dc;
    }
    rv[_nmeshcellsX-1] = h03;

    return rv;
}

void vfps::PhaseSpace::swap(vfps::PhaseSpace& other) noexcept
{
    std::swap(_data, other._data);
}

#if INOVESA_USE_OPENCL == 1
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


void vfps::swap(vfps::PhaseSpace &first, vfps::PhaseSpace &second) noexcept
{
    first.swap(second);
}
