/******************************************************************************
 * Inovesa - Inovesa Numerical Optimized Vlasov-Equation Solver Application   *
 * Copyright (c) 2014-2018: Patrik Schönfeldt                                 *
 * Copyright (c) 2014-2018: Karlsruhe Institute of Technology                 *
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

#if INOVESA_USE_OPENGL == 1

#include "IO/GUI/GUIElement.hpp"

#include <exception>

class ShaderException : public std::exception
{
public:
    ShaderException(const std::string msg) :
        _msg(msg)
                {}

        const char* what() const noexcept
        { return _msg.c_str(); }

private:
    std::string _msg;
};

vfps::GUIElement::GUIElement(bool buffer_shared)
 : buffer_shared(buffer_shared)
{
}

vfps::GUIElement::~GUIElement() noexcept
{
    glDeleteProgram(programID);
}

void vfps::GUIElement::compileShaders()
{
    // Create the shaders
    GLuint VertexShaderID = glCreateShader(GL_VERTEX_SHADER);
    GLuint FragmentShaderID = glCreateShader(GL_FRAGMENT_SHADER);

	GLint Result = GL_FALSE;
	GLint InfoLogLength;

        // Compile Vertex Shader
    char const * VertexSourcePointer = _vertexshadercode.c_str();
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
    char const * FragmentSourcePointer = _fragmentshadercode.c_str();
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

#endif // INOVESA_USE_OPENGL
