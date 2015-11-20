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

#include "IO/GUI/GUIElement.hpp"

#include <exception>

class ShaderException : public std::exception
{
public:
    ShaderException(std::string fname) :
		_fname(fname)
		{}

	const char* what() const noexcept
		{ return ("Cannot open '"+_fname+"'.").c_str(); }

private:
	std::string _fname;
};

vfps::GUIElement::GUIElement()
{
}

vfps::GUIElement::~GUIElement()
{
    glDeleteProgram(programID);
}

void vfps::GUIElement::loadShaders(std::string vertex_file_path,
                                   std::string fragment_file_path)
{
	// Create the shaders
	GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
	GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

	// Read the Vertex Shader code from the file
	std::string VertexShaderCode;
	std::ifstream VertexShaderStream(vertex_file_path, std::ios::in);
	if(VertexShaderStream.is_open()){
		std::string Line = "";
		while(getline(VertexShaderStream, Line))
			VertexShaderCode += "\n" + Line;
		VertexShaderStream.close();
	}else{
        throw ShaderException(vertex_file_path);
	}

	// Read the Fragment Shader code from the file
	std::string FragmentShaderCode;
	std::ifstream FragmentShaderStream(fragment_file_path, std::ios::in);
	if(FragmentShaderStream.is_open()){
		std::string Line = "";
		while(getline(FragmentShaderStream, Line))
			FragmentShaderCode += "\n" + Line;
		FragmentShaderStream.close();
	}

	GLint Result = GL_FALSE;
	GLint InfoLogLength;

	// Compile Vertex Shader
	char const * VertexSourcePointer = VertexShaderCode.c_str();
    glShaderSource(VertexShaderID, 1, &VertexSourcePointer , nullptr);
	glCompileShader(VertexShaderID);

	// Check Vertex Shader
	glGetShaderiv(VertexShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(VertexShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 1 ){
        char* VertexShaderErrorMessage = new char[InfoLogLength+1];
        VertexShaderErrorMessage[InfoLogLength] = '\0';
        glGetShaderInfoLog(VertexShaderID, InfoLogLength, nullptr, VertexShaderErrorMessage);
        throw ShaderException(VertexShaderErrorMessage);
	}

	// Compile Fragment Shader
	char const * FragmentSourcePointer = FragmentShaderCode.c_str();
    glShaderSource(FragmentShaderID, 1, &FragmentSourcePointer , nullptr);
	glCompileShader(FragmentShaderID);

	// Check Fragment Shader
	glGetShaderiv(FragmentShaderID, GL_COMPILE_STATUS, &Result);
	glGetShaderiv(FragmentShaderID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 1 ){
        char* FragmentShaderErrorMessage = new char[InfoLogLength+1];
        FragmentShaderErrorMessage[InfoLogLength] = '\0';
        glGetShaderInfoLog(FragmentShaderID, InfoLogLength, nullptr, FragmentShaderErrorMessage);
        throw ShaderException(FragmentShaderErrorMessage);
	}

	// Link the program
    programID = glCreateProgram();
    glAttachShader(programID, VertexShaderID);
    glAttachShader(programID, FragmentShaderID);
    glLinkProgram(programID);

	// Check the program
    glGetProgramiv(programID, GL_LINK_STATUS, &Result);
    glGetProgramiv(programID, GL_INFO_LOG_LENGTH, &InfoLogLength);
	if ( InfoLogLength > 1 ){
        char* ProgramErrorMessage = new char[InfoLogLength+1];
        ProgramErrorMessage[InfoLogLength] = '\0';
        glGetProgramInfoLog(programID, InfoLogLength, nullptr, &ProgramErrorMessage[0]);
        throw ShaderException(ProgramErrorMessage);
	}

	glDeleteShader(VertexShaderID);
    glDeleteShader(FragmentShaderID);
}

uint_fast8_t vfps::GUIElement::glversion;
