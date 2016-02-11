/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Algorithms   *
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

#include "IO/GUI/Plot3DColormap.hpp"

vfps::Plot3DColormap::Plot3DColormap(meshdata_t maxvalue) :
    maxValue(maxvalue)
{
    // Create GLSL code
    switch (glversion) {
    case 2:
        _fragmentshadercode = R"(
            #version 120

            // Interpolated values from the vertex shaders
            varying vec2 UV3DCM;

            // Values that stay constant for the whole mesh.
            uniform sampler2D textureSampler3DCM;

            void main(){

                // Output color = color of the texture at the specified UV
                gl_FragColor =  texture2D( textureSampler3DCM, UV3DCM );
            }
            )";
        _vertexshadercode = R"(
            #version 120

            // Input vertex data, different for all executions of this shader.
            attribute vec2 position3DCM;
            attribute vec2 vertexUV3DCM;

            // Output data ; will be interpolated for each fragment.
            varying vec2 UV3DCM;

            void main(){

                gl_Position =  vec4(position3DCM,0,1);

                // UV of the vertex. No special space for this one.
                UV3DCM = vertexUV3DCM;
            }
            )";
        break;
    case 3:
    default:
        _fragmentshadercode = R"(
            #version 330 core

            // Interpolated values from the vertex shaders
            in vec2 UV3DCM;

            // Ouput data
            out vec3 color3DCM;

            // Values that stay constant for the whole mesh.
            uniform sampler2D textureSampler3DCM;

            void main(){
                // Output color = color of the texture at the specified UV
                color3DCM = texture2D( textureSampler3DCM, UV3DCM ).rgb;
            }
            )";
        _vertexshadercode = R"(
            #version 330 core

            // Input vertex data, different for all executions of this shader.
            layout(location = 0) in vec2 position3DCM;
            layout(location = 1) in vec2 vertexUV3DCM;

            // Output data ; will be interpolated for each fragment.
            out vec2 UV3DCM;

            void main(){

                // Output position of the vertex, in clip space : MVP * position
                gl_Position = vec4(position3DCM,0,1);

                // UV of the vertex. No special space for this one.
                UV3DCM = vertexUV3DCM;
            }
            )";
        break;
    }
    compileShaders();
    position = glGetAttribLocation(programID, "position3DCM");
    vertexUV = glGetAttribLocation(programID, "vertexUV3DCM");
    textureSampler = glGetUniformLocation(programID, "textureSampler3DCM");

    static const GLfloat g_vertex_buffer_data[] = {
        -1.0f,-0.5f,
         0.5f,-0.5f,
        -1.0f, 1.0f,
         0.5f, 1.0f,
        -1.0f, 1.0f,
         0.5f,-0.5f
    };

	/* Two UV coordinatesfor each vertex.
	 * Cordinates are switched because Mesh2D has data arranged as [x][y],
	 * which is not the way OpenGL would expect it.
	 */
    static const GLfloat g_uv_buffer_data[] = {
        0.0f, 0.0f,
        0.0f, 1.0f,
        1.0f, 0.0f,
        1.0f, 1.0f,
        1.0f, 0.0f,
        0.0f, 1.0f
    };
        glGenBuffers(1, &vertexbuffer);
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(g_vertex_buffer_data),
                                 g_vertex_buffer_data, GL_STATIC_DRAW);

	glGenBuffers(1, &uvbuffer);
	glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
	glBufferData(GL_ARRAY_BUFFER, sizeof(g_uv_buffer_data),
				 g_uv_buffer_data, GL_STATIC_DRAW);
}

vfps::Plot3DColormap::~Plot3DColormap()
{
    glDeleteBuffers(1,&uvbuffer);
    glDeleteBuffers(1,&vertexbuffer);
    glDeleteVertexArrays(1,&vertexUV);
    glDeleteVertexArrays(1,&position);
    glDeleteTextures(1, &textureID);
    glDeleteTextures(1, &textureSampler);
}


void vfps::Plot3DColormap::createTexture(vfps::PhaseSpace* mesh)
{
    glGenTextures (1, &textureID);
        glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
    glBindTexture(GL_TEXTURE_2D,textureID);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);

	size_t npixels = mesh->nMeshCells();
	float* data = new float[npixels];
	vfps::meshdata_t* meshdata = mesh->getData();
	for (vfps::meshindex_t i=0; i<npixels; i++) {
		data[i] = meshdata[i]/maxValue;
	}
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F,
                 mesh->nMeshCells(0), mesh->nMeshCells(1),
                 0, GL_RED, GL_FLOAT, data);
        delete [] data;

	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,  GL_CLAMP_TO_BORDER );
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,  GL_CLAMP_TO_BORDER );
    glBindTexture(GL_TEXTURE_2D, 0);
}

void vfps::Plot3DColormap::delTexture()
{
    glDeleteTextures(1, &textureID);
}

void vfps::Plot3DColormap::draw()
{
	// Use our shader
	glUseProgram(programID);

	// Bind our texture in Texture Unit 0
	glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, textureID);
    // Set "textureSampler3DCM" sampler to user Texture Unit 0
    glUniform1i(textureSampler, 0);

    glEnableVertexAttribArray(position);
        glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(position,2,GL_FLOAT,GL_FALSE,0,nullptr);

    glEnableVertexAttribArray(vertexUV);
    glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
    glVertexAttribPointer(vertexUV,2,GL_FLOAT,GL_FALSE,0,nullptr);


    glDrawArrays(GL_TRIANGLES, 0, 2*3);// 2*3 indices starting at 0

    glDisableVertexAttribArray(position);
    glDisableVertexAttribArray(vertexUV);
}

#endif // INOVESA_USE_GUI
