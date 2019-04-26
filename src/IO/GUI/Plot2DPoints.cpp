// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#if INOVESA_USE_OPENGL == 1

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
    _data.resize(2*points.size());

    for (size_t n=0; n<_npoints; n++) {
        _data[2*n  ] = points[n].x*1.5f/_nmeshcellsx-1.0f;
        _data[2*n+1] = points[n].y*1.5f/_nmeshcellsy-0.5f;
    }


    // The following commands will talk about our 'vertexbuffer' buffer
    glBindBuffer(GL_ARRAY_BUFFER, databuffer);

    // Give our vertices to OpenGL.
    glBufferData(GL_ARRAY_BUFFER, 2*_npoints*sizeof(float),
                 _data.data(), GL_STATIC_DRAW);
}

#endif // INOVESA_USE_OPENGL
