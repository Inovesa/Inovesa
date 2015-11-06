#include "IO/GUI/Plot2DLine.hpp"

vfps::Plot2DLine::Plot2DLine() :
    _npoints(3),
    _line(new float[3*_npoints])
{
    // Create and compile our GLSL program from the shaders
    switch (glversion) {
    case 2:
        programID = LoadShaders("gl/2DL_gl2.vertexshader",
                                "gl/2DL_gl2.fragmentshader");
        break;
    case 3:
    default:
        programID = LoadShaders("gl/2DL_gl3.vertexshader",
                                "gl/2DL_gl3.fragmentshader");
        break;
    }
    VertexArrayID = glGetAttribLocation(programID, "Position2DL");

    createLine(nullptr);

    // Generate 1 buffer, put the resulting identifier in vertexbuffer
    glGenBuffers(1, &vertexbuffer);

    // The following commands will talk about our 'vertexbuffer' buffer
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);

    // Give our vertices to OpenGL.
    glBufferData(GL_ARRAY_BUFFER, 3*_npoints*sizeof(float),
                 _line, GL_STATIC_DRAW);
}

vfps::Plot2DLine::~Plot2DLine()
{
    delete [] _line;
    glDeleteProgram(programID);
    glDeleteBuffers(1,&vertexbuffer);
}

void vfps::Plot2DLine::createLine(vfps::PhaseSpace* mesh)
{
    float step = 0.1f;
    for (size_t n=0; n<3*_npoints; n+=3) {
        _line[n  ] = -0.9f+n*step;
        _line[n+1] = -0.8f;
        _line[n+2] =  0.0f;
    }
}

void vfps::Plot2DLine::draw()
{
    glUseProgram(programID);

    // 1rst attribute buffer : vertices
    glEnableVertexAttribArray(VertexArrayID);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(
       0,                  // attribute 0. No particular reason for 0,
                           // but must match the layout in the shader.
       _npoints,           // size
       GL_FLOAT,           // type
       GL_FALSE,           // normalized?
       0,                  // stride
       (void*)0            // array buffer offset
    );

    // Draw the triangle !GL_LINE_STRIP
    glDrawArrays(GL_LINE_STRIP, 0, _npoints); // Starting from vertex 0; 3 vertices total -> 1 triangle

    glDisableVertexAttribArray(VertexArrayID);
}
