/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
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

#ifdef INOVESA_USE_GUI

#include "IO/GUI/Plot2DLine.hpp"

vfps::Plot2DLine::Plot2DLine(std::array<float,3> rgb) :
    _npoints(0)
{
    std::stringstream color;
    color << rgb[0] << ',' << rgb[1] << ',' << rgb[2];
    // Create GLSL code
    switch (glversion) {
    case 2:
        _fragmentshadercode = R"(
                #version 120
                void main(){
                    gl_FragColor = vec4(
        )";
        _fragmentshadercode += color.str() + ",0);}";
        _vertexshadercode = R"(
            #version 120

            // Input vertex data, different for all executions of this shader.
            attribute vec2 position2DL;

            void main(){
                gl_Position =  vec4(position2DL,0,1);
            }
            )";
        break;
    case 3:
    default:
        _fragmentshadercode = R"(
            #version 330 core
            out vec3 color2DL;

            void main(){
                color2DL = vec3(
            )";
            _fragmentshadercode += color.str() + ");}";
        _vertexshadercode = R"(
            #version 330 core

            // Input vertex data, different for all executions of this shader.
            layout(location = 0) in vec2 position2DL;

            void main(){
                gl_Position =  vec4(position2DL,0,1);
            }
            )";
        break;
    }
    compileShaders();
    position = glGetAttribLocation(programID, "position2DL");

    // Generate 1 buffer, put the resulting identifier in vertexbuffer
    glGenBuffers(1, &vertexbuffer);
}

vfps::Plot2DLine::~Plot2DLine()
{
    glDeleteBuffers(1, &vertexbuffer);
    glDeleteVertexArrays(1, &position);
}

void vfps::Plot2DLine::updateLine(const size_t npoints,
                                  const float* points,
                                  const bool vertical)
{
    _npoints=npoints;
    _line.resize(_npoints*2);
    float step = 1.5f/_npoints;
    float max =0.0f;
    for (size_t n=0; n<_npoints; n++) {
        max = std::max(max,points[n]);
    }
    if (vertical) {
        max*=2.0f;
        for (size_t n=0; n<_npoints; n++) {
            _line[2*n  ] = 0.5f+points[n]/max;
            _line[2*n+1] = -0.5f+n*step;
        }
    } else {
        max*=4.0f;
        for (size_t n=0; n<_npoints; n++) {
            _line[2*n  ] = -1.0f+n*step;
            _line[2*n+1] = -0.8f+points[n]/max;
        }
    }

    // The following commands will talk about our 'vertexbuffer' buffer
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);

    // Give our vertices to OpenGL.
    glBufferData(GL_ARRAY_BUFFER, 2*_npoints*sizeof(float),
                 _line.data(), GL_STATIC_DRAW);
}

void vfps::Plot2DLine::draw()
{
    glUseProgram(programID);

    glEnableVertexAttribArray(position);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(position,2,GL_FLOAT,GL_FALSE,0,nullptr);

    // Draw the line
    glDrawArrays(GL_LINE_STRIP, 0, _npoints);

    glDisableVertexAttribArray(position);
}

#endif // INOVESA_USE_GUI
