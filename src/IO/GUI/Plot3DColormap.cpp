#include "IO/GUI/Plot3DColormap.hpp"

vfps::Plot3DColormap::Plot3DColormap()
{
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	switch (glversion) {
	case 2:
		programID = LoadShaders("gl/gl2.vertexshader","gl/gl2.fragmentshader");
		break;
	case 3:
	default:
		programID = LoadShaders("gl/gl3.vertexshader","gl/gl3.fragmentshader");
		break;
	}

	// Get a handle for our "MVP" uniform
	MatrixID = glGetUniformLocation(programID, "MVP");

	// Projection matrix : 45° Field of View, 1:1 ratio, display range : 0.1 unit <-> 100 units
	glm::mat4 Projection = glm::perspective(45.0f, 1.0f, 0.1f, 100.0f);
	// Camera matrix
	glm::mat4 View	   = glm::lookAt(
								glm::vec3(0.5,0.5,1.0),
								glm::vec3(0.5,0.5,0),
								glm::vec3(0,1,0)
						   );
	// Model matrix : an identity matrix (model will be at the origin)
	glm::mat4 Model	  = glm::mat4(1.0f);
	// Our ModelViewProjection : multiplication of our 3 matrices
	MVP  = Projection * View * Model;

	static const GLfloat g_vertex_buffer_data[] = {
		0.0f, 0.0f, 0.0f,
		1.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
		1.0f, 1.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
		1.0f, 0.0f, 0.0f
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

	// Get a handle for our "myTextureSampler" uniform
	TextureID = glGetUniformLocation(programID, "myTextureSampler");
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

	// Send our transformation to the currently bound shader,
	// in the "MVP" uniform
	glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);

	// Bind our texture in Texture Unit 0
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, Texture);
	// Set our "myTextureSampler" sampler to user Texture Unit 0
	glUniform1i(TextureID, 0);

	// 1rst attribute buffer : vertices
	glEnableVertexAttribArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
	glVertexAttribPointer(
		0,				  // attribute. No particular reason for 0, but must match the layout in the shader.
		3,				  // size
		GL_FLOAT,		   // type
		GL_FALSE,		   // normalized?
		0,				  // stride
		(void*)0			// array buffer offset
	);

	// 2nd attribute buffer : UVs
	glEnableVertexAttribArray(1);
	glBindBuffer(GL_ARRAY_BUFFER, uvbuffer);
	glVertexAttribPointer(
		1,								// attribute. No particular reason for 1, but must match the layout in the shader.
		2,								// size : U+V => 2
		GL_FLOAT,						 // type
		GL_FALSE,						 // normalized?
		0,								// stride
		(void*)0						  // array buffer offset
	);
	// Draw the triangle !
	glDrawArrays(GL_TRIANGLES, 0, 2*3); // 12*3 indices starting at 0 -> 12 triangles

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
}
