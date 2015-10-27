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

vfps::Display::Display()
	#if GLFW_VERSION_MAJOR == 3
	:
	window(nullptr)
	#endif
{
	glfwInit();

	#if GLFW_VERSION_MAJOR == 3 // GLFW3
	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	#else // GLFW2
	glfwOpenWindowHint(GLFW_FSAA_SAMPLES, 4);
	glfwOpenWindowHint(GLFW_OPENGL_VERSION_MAJOR, 3);
	glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 3);
	glfwOpenWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	#endif // GLFW2


	// Open a window and create its OpenGL context
	#if GLFW_VERSION_MAJOR < 3
	glfwOpenWindow( 512, 512,6,5,6,0,0,0, GLFW_WINDOW);
	glfwSetWindowTitle("Inovesa");
	#else // GLFW3
	window = glfwCreateWindow( 512, 512, "Inovesa", NULL, NULL);
	if( window == nullptr ) {
		GUIElement::glversion = 2;

		glfwTerminate();
		printText("Failed to open OpenGl 3 window, will try OpenGL 2.");
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
		GUIElement::glversion = 3;
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
}

vfps::Display::~Display()
{
	// Close OpenGL window and terminate GLFW
	glfwTerminate();
}

void vfps::Display::addElement(GUIElement* newitem)
{
	_item.push_back(newitem);
}

void vfps::Display::draw() {
	// Clear the screen
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// draw GUIElements
	for (GUIElement* i : _item) {
		i->draw();
	}

	// Swap buffers
	#if GLFW_VERSION_MAJOR == 3
	glfwSwapBuffers(window);
	#else
	glfwSwapBuffers();
	#endif
	glfwPollEvents();
}

#endif // INOVESA_USE_GUI

void vfps::Display::printText(std::string txt, bool signOfLife)
{
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    if (!signOfLife ||
        std::chrono::duration<double>(now-_lastmessage).count() > 1.0) {
        std::cout.setf( std::ios::fixed, std:: ios::floatfield );
        std::cout.precision(3);
        std::cout << "[ " << std::setw(9)
                  << std::chrono::duration<double>(now-start_time).count()
                  << " ]: "
                  << txt
                  << std::endl;
        _lastmessage = now;
    }
}

void vfps::Display::takeElement(vfps::GUIElement* item)
{
	for (size_t i=0; i< _item.size(); i++) {
		if (_item[i] == item) {
			_item.erase(_item.begin()+i);
			i--;
		}
	}
}

std::chrono::system_clock::time_point vfps::Display::start_time;

std::chrono::system_clock::time_point vfps::Display::_lastmessage;

