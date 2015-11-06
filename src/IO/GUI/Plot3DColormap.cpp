#include "IO/GUI/Plot3DColormap.hpp"

vfps::Plot3DColormap::Plot3DColormap()
{
	// Create and compile our GLSL program from the shaders
	switch (glversion) {
	case 2:
        loadShaders("gl/3DCM_gl2.vertexshader","gl/3DCM_gl2.fragmentshader");
		break;
	case 3:
	default:
        loadShaders("gl/3DCM_gl3.vertexshader","gl/3DCM_gl3.fragmentshader");
		break;
	}
    position = glGetAttribLocation(programID, "position3DCM");
    vertexUV = glGetAttribLocation(programID, "vertexUV3DCM");
    TextureID = glGetUniformLocation(programID, "textureSampler3DCM");

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
    delTexture();
    glDeleteBuffers(1,&uvbuffer);
    glDeleteBuffers(1,&vertexbuffer);
    glDeleteVertexArrays(1,&vertexUV);
    glDeleteVertexArrays(1,&position);
}


void vfps::Plot3DColormap::createTexture(vfps::PhaseSpace* mesh)
{
	glGenTextures (1, &Texture);
	glTexEnvi( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE );
	glBindTexture(GL_TEXTURE_2D,Texture);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);

	size_t npixels = mesh->nMeshCells();
	float* data = new float[npixels];
	vfps::meshdata_t* meshdata = mesh->getData();
	for (vfps::meshindex_t i=0; i<npixels; i++) {
		data[i] = meshdata[i]/vfps::meshdata_t(1);
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
    glDeleteTextures(1, &TextureID);
    glDeleteTextures(1, &Texture);
}

void vfps::Plot3DColormap::draw()
{
	// Use our shader
	glUseProgram(programID);

	// Bind our texture in Texture Unit 0
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, Texture);
    // Set "textureSampler3DCM" sampler to user Texture Unit 0
    glUniform1i(TextureID, 0);

    glEnableVertexAttribArray(position);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
    glVertexAttribPointer(position,2,GL_FLOAT,GL_FALSE,0,nullptr);


    glEnableVertexAttribArray(vertexUV);
	glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
    glVertexAttribPointer(1,2,GL_FLOAT,GL_FALSE,0,nullptr);

    glDrawArrays(GL_TRIANGLES, 0, 2*3);// 2*3 indices starting at 0

    glDisableVertexAttribArray(position);
    glDisableVertexAttribArray(vertexUV);
}
