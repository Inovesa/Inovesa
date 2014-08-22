#ifndef DISPLAY_HPP
#define DISPLAY_HPP

#include <array>
#include <fstream>
#include <iostream>
#include <type_traits>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#ifdef GLFW3
#include <GLFW/glfw3.h>
#else // GLFW3
#include <GL/glfw.h>
#endif // GLFW3

// Include GLM
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "Mesh2D.hpp"

class Display
{
public:
	Display();

	~Display();

	void draw(vfps::Mesh2D<meshdata_t>* mesh);

private:
	GLuint LoadShaders(const char* vertex_file_path,const char* fragment_file_path);

	#ifdef GLFW3
	GLFWwindow* window;
	#endif
	bool gl2fallback;

	GLuint vertexbuffer;
	GLuint uvbuffer;
	GLuint programID;
	GLuint VertexArrayID;
	GLuint MatrixID;
	glm::mat4 MVP;
};

#endif // DISPLAY_HPP
