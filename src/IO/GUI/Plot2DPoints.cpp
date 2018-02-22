/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2018: Patrik Schönfeldt                                      *
 * Copyright (c) 2018: Karlsruhe Institute of Technology                      *
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

#include "IO/GUI/Plot2DPoints.hpp"


vfps::Plot2DPoints::Plot2DPoints(std::array<float,3> rgb,
                                 uint32_t nmeshcellsx,
                                 uint32_t nmeshcellsy)
 : Plot2DPrimitive(rgb,0,GL_POINTS)
 , _nmeshcellsx(nmeshcellsx)
 , _nmeshcellsy(nmeshcellsy)
{
}

void vfps::Plot2DPoints::update(const std::vector<PhaseSpace::Position> &points)
{
    _npoints=points.size();
    _points.resize(2*points.size());

    for (size_t n=0; n<_npoints; n++) {
        _points[2*n  ] = points[n].x*1.5f/_nmeshcellsx-1.0f;
        _points[2*n+1] = points[n].y*1.5f/_nmeshcellsy-0.5f;
    }


    // The following commands will talk about our 'vertexbuffer' buffer
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);

    // Give our vertices to OpenGL.
    glBufferData(GL_ARRAY_BUFFER, 2*_npoints*sizeof(float),
                 _points.data(), GL_STATIC_DRAW);
}

#endif // INOVESA_USE_OPENGL
