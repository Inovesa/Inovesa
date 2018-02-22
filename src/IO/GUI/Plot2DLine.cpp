/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2018: Patrik Schönfeldt                                 *
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

#ifdef INOVESA_USE_OPENGL

#include "IO/GUI/Plot2DLine.hpp"

vfps::Plot2DLine::Plot2DLine(std::array<float,3> rgb)
  : Plot2DPrimitive(rgb,0,GL_LINE_STRIP)
  , _max(0)
{
}

void vfps::Plot2DLine::update(const size_t npoints,
                              const float* points,
                              const bool vertical)
{
    _npoints=npoints;
    _points.resize(npoints*2);
    float step = 1.5f/(_npoints-1);
    for (size_t n=0; n<_npoints; n++) {
        _max = std::max(_max,std::abs(points[n]));
    }
    if (vertical) {
        const float max = 2*_max;
        for (size_t n=0; n<npoints; n++) {
            _points[2*n  ] = 0.5f+points[n]/max;
            _points[2*n+1] = -0.5f+n*step;
        }
    } else {
        const float max = 4*_max;
        for (size_t n=0; n<_npoints; n++) {
            _points[2*n  ] = -1.0f+n*step;
            _points[2*n+1] = -0.8f+points[n]/max;
        }
    }

    // The following commands will talk about our 'vertexbuffer' buffer
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);

    // Give our vertices to OpenGL.
    glBufferData(GL_ARRAY_BUFFER, 2*_npoints*sizeof(float),
                 _points.data(), GL_STATIC_DRAW);
}

void vfps::Plot2DLine::update(const size_t npoints,
                              const double* points,
                              const bool vertical)
{
    std::vector<float> fltpoints(npoints);
    for (size_t i=0; i<npoints; i++) {
        fltpoints[i] = points[i];
    }
    update(npoints,fltpoints.data(),vertical);
}

#endif // INOVESA_USE_OPENGL
