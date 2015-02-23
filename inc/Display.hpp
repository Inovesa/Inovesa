/******************************************************************************/
/* Inovesa - Inovesa Numerical Optimized Vlesov-Equation Solver Application   */
/* Copyright (c) 2007-2009: Peter Schregle (Fixed Point Math Library)         */
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

#include <array>
#include <fstream>
#include <iostream>
#include <type_traits>

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#define GLFW3
#ifdef GLFW3
#include <GLFW/glfw3.h>
#else // GLFW3
#include <GL/glfw.h>
#endif // GLFW3

// Include GLM
#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include "PhaseSpace.hpp"

class Display
{
public:
	Display();

	~Display();

	void createTexture(vfps::PhaseSpace* mesh);

	void delTexture();

	void draw();

	GLuint getTexture() const
		{ return Texture; }

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
	GLuint Texture;
	GLuint TextureID;
	glm::mat4 MVP;
};

#endif // DISPLAY_HPP
