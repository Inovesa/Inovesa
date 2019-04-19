// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#if INOVESA_USE_OPENGL == 1

#include "IO/GUI/Plot2DPrimitive.hpp"

vfps::Plot2DPrimitive::Plot2DPrimitive(std::array<float,3> rgb,
                                       size_t npoints,
                                       GLenum primitivetype)
  : GUIElement(false)
  , _npoints(npoints)
  , _primitivetype(primitivetype)
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
                gl_Position =  vec4(position2DL,-0.1,1);
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
                gl_Position =  vec4(position2DL,-0.1,1);
            }
            )";
        break;
    }
    compileShaders();
    dataID = glGetAttribLocation(programID, "position2DL");

    // Generate 1 buffer, put the resulting identifier in databuffer
    glGenBuffers(1, &databuffer);
}

vfps::Plot2DPrimitive::~Plot2DPrimitive() noexcept
{
    glDeleteBuffers(1, &databuffer);
    glDeleteVertexArrays(1, &dataID);
    glDeleteProgram(programID);
}

void vfps::Plot2DPrimitive::draw()
{
    glUseProgram(programID);

    glEnableVertexAttribArray(dataID);
    glBindBuffer(GL_ARRAY_BUFFER, databuffer);
    glVertexAttribPointer(dataID,2,GL_FLOAT,GL_FALSE,0,nullptr);

    // Draw the line
    glDrawArrays(_primitivetype, 0, _npoints);

    glDisableVertexAttribArray(dataID);
}

#endif // INOVESA_USE_OPENGL
