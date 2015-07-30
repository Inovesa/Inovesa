/******************************************************************************/
/* Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   */
/* Copyright (c) 2014-2015: Patrik Schönfeldt                                 */
/*                                                                            */
/* This file is part of Inovesa.                                              */
/* Inovesa is free software: you can redistribute it and/or modify            */
/* it under the terms of the GNU General Public License as published by       */
/* the Free Software Foundation, either version 3 of the License, or          */
/* (at your option) any later version.                                        */
/*                                                                            */
/* Inovesa is distributed in the hope that it will be useful,                 */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with Inovesa.  If not, see <http://www.gnu.org/licenses/>.           */
/******************************************************************************/

#include "IO/Display.hpp"

#ifdef INOVESA_USE_GUI

vfps::Display::Display() :
	#if GLFW_VERSION_MAJOR == 3
	window(nullptr)
	#endif
{
	glfwInit();

	#if GLFW_VERSION_MAJOR == 3
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	#else // GLFW3
	glfwOpenWindowHint(GLFW_FSAA_SAMPLES, 4);
	glfwOpenWindowHint(GLFW_OPENGL_VERSION_MAJOR, 3);
	glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 3);
	glfwOpenWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	#endif // GLFW3


	// Open a window and create its OpenGL context
	#if GLFW_VERSION_MAJOR < 3
	glfwOpenWindow( 512, 512,6,5,6,0,0,0, GLFW_WINDOW);
	glfwSetWindowTitle("Inovesa");
	#else // GLFW3
	window = glfwCreateWindow( 512, 512, "Inovesa", NULL, NULL);
	if( window == nullptr ) {
		glfwTerminate();
		printText("Failed to open OpenGl 3 window, will try OpenGL 2.");
		gl2fallback = true;
		glfwInit();
		glfwWindowHint(GLFW_SAMPLES, 4);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
		glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
		window = glfwCreateWindow( 512, 512, "Phace Space View", NULL, NULL);

		if( window == nullptr ) {
			std::cerr << "Failed to initialize GLFW." << std::endl;
			glfwTerminate();
		}
	} else {
		gl2fallback = false;
	}
	glfwMakeContextCurrent(window);
	#endif // GLFW3

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		std::cerr << "Failed to initialize GLEW" << std::endl;
	}

	// Ensure we can capture the escape key being pressed below
	#if GLFW_VERSION_MAJOR == 3
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);
	#else
	glfwEnable(GLFW_STICKY_KEYS);
	#endif

	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);

	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

	// Create and compile our GLSL program from the shaders
	if (!gl2fallback) {
		programID = LoadShaders("gl/gl3.vertexshader","gl/gl3.fragmentshader");
	} else {
		programID = LoadShaders("gl/gl2.vertexshader","gl/gl2.fragmentshader");
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

vfps::Display::~Display()
{
	delTexture();

	// Cleanup VBO and shader
	glDeleteBuffers(1, &vertexbuffer);
	glDeleteBuffers(1, &uvbuffer);
	glDeleteProgram(programID);
	glDeleteVertexArrays(1, &VertexArrayID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();
}

void vfps::Display::draw() {
	// Clear the screen
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

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

	// Swap buffers
	#if GLFW_VERSION_MAJOR == 3
	glfwSwapBuffers(window);
	#else
	glfwSwapBuffers();
	#endif
	glfwPollEvents();
}

#endif // INOVESA_USE_GUI

void vfps::Display::printText(std::string txt)
{
	std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
	std::cout.setf( std::ios::fixed, std:: ios::floatfield );
	std::cout.precision(3);
	std::cout
				<< "[ " << std::setw(9)
				<< std::chrono::duration<double>(now-start_time).count()
				<< " ]: "
				<< txt
				<< std::endl;
}

std::chrono::system_clock::time_point vfps::Display::start_time;

bool vfps::Display::gl2fallback;

