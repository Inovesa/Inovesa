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

#include "IO/GUI/Plot1DLine.hpp"

vfps::Plot1DLine::Plot1DLine(std::array<float,3> rgb
                            , size_t npoints
                            , orientation orientation)
  : _npoints(npoints)
  , _max(0)
  , _orientation(orientation)
{
    _data.resize(npoints);
    _position.resize(npoints);

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
            attribute float posiX;
            attribute float posiY;

            void main(){
                gl_Position =  vec4(posiX,posiY,-0.1,1);
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
            layout(location = 0) in vec posiX;
            layout(location = 1) in vec posiY;

            void main(){
                gl_Position =  vec4(posiX,posiY,-0.1,1);
            }
            )";
        break;
    }
    compileShaders();

    switch (_orientation) {
    case orientation::horizontal:
        posiID = glGetAttribLocation(programID, "posiX");
        dataID = glGetAttribLocation(programID, "posiY");
        break;
    case orientation::vertical:
        posiID = glGetAttribLocation(programID, "posiY");
        dataID = glGetAttribLocation(programID, "posiX");
        break;
    }

    glGenBuffers(1, &databuffer);
    glGenBuffers(1, &posibuffer);

    createPositionBuffer();
}

vfps::Plot1DLine::~Plot1DLine() noexcept
{
    glDeleteBuffers(1, &databuffer);
    glDeleteVertexArrays(1, &dataID);
    glDeleteBuffers(1, &posibuffer);
    glDeleteVertexArrays(1, &posiID);
    glDeleteProgram(programID);
}

void vfps::Plot1DLine::draw()
{
    glUseProgram(programID);

    glEnableVertexAttribArray(posiID);
    glBindBuffer(GL_ARRAY_BUFFER, posibuffer);
    glVertexAttribPointer(posiID,1,GL_FLOAT,GL_FALSE,0,nullptr);

    glEnableVertexAttribArray(dataID);
    glBindBuffer(GL_ARRAY_BUFFER, databuffer);
    glVertexAttribPointer(dataID,1,GL_FLOAT,GL_FALSE,0,nullptr);

    // Draw the line
    glDrawArrays(GL_LINE_STRIP, 0, _npoints);

    glDisableVertexAttribArray(posiID);
    glDisableVertexAttribArray(dataID);
}

void vfps::Plot1DLine::update(const size_t npoints,
                              const float* points)
{

    for (size_t n=0; n<_npoints; n++) {
        _max = std::max(_max,std::abs(points[n]));
    }
    switch (_orientation) {
    case orientation::horizontal:
        {
        const float max = 4*_max;
        for (size_t n=0; n<_npoints; n++) {
            _data[n] = -0.8f+points[n]/max;
        }
        }
        break;
    case orientation::vertical:
        {
        const float max = 2*_max;
        for (size_t n=0; n<npoints; n++) {
            _data[n] = 0.5f+points[n]/max;
        }
        break;
        }
    }

    glBindBuffer(GL_ARRAY_BUFFER, databuffer);
    glBufferData( GL_ARRAY_BUFFER, _npoints*sizeof(float)
                , _data.data(), GL_STATIC_DRAW);
}

void vfps::Plot1DLine::update(const size_t npoints,
                              const double* points)
{
    std::vector<float> fltpoints(npoints);
    for (size_t i=0; i<npoints; i++) {
        fltpoints[i] = points[i];
    }
    update(npoints,fltpoints.data());
}

void vfps::Plot1DLine::createPositionBuffer()
{
    float step = 1.5f/(_npoints-1);
    switch (_orientation) {
    case orientation::horizontal:
        for (size_t n=0; n<_npoints; n++) {
            _position[n] = -1.0f+n*step;
        }
        break;
    case orientation::vertical:
        for (size_t n=0; n<_npoints; n++) {
            _position[n] = -0.5f+n*step;;
        }
        break;
    }

    glBindBuffer(GL_ARRAY_BUFFER, posibuffer);
    glBufferData( GL_ARRAY_BUFFER, _npoints*sizeof(float)
                , _position.data(), GL_STATIC_DRAW);
}

#endif // INOVESA_USE_OPENGL
