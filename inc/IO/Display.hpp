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
#include <iomanip>
#include <iostream>
#include <type_traits>

#ifdef INOVESA_USE_GUI

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#if GLFW_VERSION_MAJOR == 3 // GLFW3
#include <GLFW/glfw3.h>
#else // GLFW2
#include <GL/glfw.h>
#endif // GLFW2
#endif // INOVESA_USE_GUI
#include <vector>

#include "PhaseSpace.hpp"
#include "IO/GUI/GUIElement.hpp"

namespace vfps {

class Display
{
public:
	static std::chrono::system_clock::time_point start_time;

public:
	Display();

	~Display();

    #ifdef INOVESA_USE_GUI
	void addElement(GUIElement* newitem);
    #endif // INOVESA_USE_GUI

	void draw();

	static void printText(std::string txt, float silentTime=0.0f);

    #ifdef INOVESA_USE_GUI
	void takeElement(GUIElement* item);
    #endif // INOVESA_USE_GUI

private:
	#ifdef INOVESA_USE_GUI
	#if GLFW_VERSION_MAJOR == 3
	GLFWwindow* window;
    #endif

	std::vector<GUIElement*> _item;
    #endif // INOVESA_USE_GUI

    static std::chrono::system_clock::time_point _lastmessage;
};

} // namespace vfps

#endif // DISPLAY_HPP

