// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

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
