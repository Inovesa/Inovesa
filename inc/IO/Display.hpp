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

#ifndef DISPLAY_HPP
#define DISPLAY_HPP

#include "defines.hpp"

#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <type_traits>

#ifdef INOVESA_USE_GUI

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#if GLFW_VERSION_MAJOR == 3
#include <GLFW/glfw3.h>
#else
#include <GL/glfw.h>
#endif

// Include GLM
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "PhaseSpace.hpp"

#endif // INOVESA_USE_GUI
class Display
{
public:
	Display();

	~Display();

	void createTexture(vfps::PhaseSpace* mesh);

	void delTexture();

	void draw();

	static void printText(std::string txt);

#ifdef INOVESA_USE_GUI
	GLuint getTexture() const
		{ return Texture; }

private:
	GLuint LoadShaders(const char* vertex_file_path,const char* fragment_file_path);

	#if GLFW_VERSION_MAJOR == 3
	GLFWwindow* window;
	#endif
	bool gl2fallback;

	GLuint vertexbuffer;
	GLuint uvbuffer;
	GLuint programID;
	GLuint VertexArrayID;
	GLuint MatrixID;
	GLuint Texture;
	GLuint TextureID;
	glm::mat4 MVP;
#endif // INOVESA_USE_GUI
};

#endif // DISPLAY_HPP

