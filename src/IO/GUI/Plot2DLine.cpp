#include "IO/GUI/Plot2DLine.hpp"

vfps::Plot2DLine::Plot2DLine() :
    _npoints(0)
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

    // Generate 1 buffer, put the resulting identifier in vertexbuffer
    glGenBuffers(1, &vertexbuffer);
}

vfps::Plot2DLine::~Plot2DLine()
{
    glDeleteProgram(programID);
    glDeleteBuffers(1,&vertexbuffer);
}

void vfps::Plot2DLine::createLine(vfps::PhaseSpace* mesh)
{
    _npoints=mesh->nMeshCells(0);
    _line.resize(_npoints*3);
    float step = 1.5f/_npoints;
    projection_t* bp = mesh->projectionToX();
    float max =0.0f;
    for (size_t n=0; n<_npoints; n++) {
        max = std::max(max,bp[n]);
    }
    max*=3.0f;
    for (size_t n=0; n<_npoints; n++) {
        _line[3*n  ] = -1.0f+n*step;
        _line[3*n+1] = -0.9f+bp[n]/max;
        _line[3*n+2] =  0.0f;
    }

    // The following commands will talk about our 'vertexbuffer' buffer
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);

    // Give our vertices to OpenGL.
    glBufferData(GL_ARRAY_BUFFER, 3*_npoints*sizeof(float),
                 _line.data(), GL_STATIC_DRAW);
}

void vfps::Plot2DLine::draw()
{
    glUseProgram(programID);

    // 1rst attribute buffer : vertices
    glEnableVertexAttribArray(VertexArrayID);
    glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(
       VertexArrayID,      // attribute 0. No particular reason for 0,
                           // but must match the layout in the shader.
       3,                  // size
       GL_FLOAT,           // type
       GL_FALSE,           // normalized?
       0,                  // stride
       (void*)0            // array buffer offset
    );

    // Draw the triangle !GL_LINE_STRIP
    glDrawArrays(GL_LINE_STRIP, 0, _npoints); // Starting from vertex 0; 3 vertices total -> 1 triangle

    glDisableVertexAttribArray(VertexArrayID);
}
