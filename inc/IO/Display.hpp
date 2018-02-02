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

#ifndef DISPLAY_HPP
#define DISPLAY_HPP

#include "defines.hpp"

#include <array>
#include <chrono>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <type_traits>

#ifdef INOVESA_USE_OPENGL

// Include GLEW
#include <GL/glew.h>

// Include GLFW
#if GLFW_VERSION_MAJOR == 3 // GLFW3
#include <GLFW/glfw3.h>
#else // GLFW2
#include <GL/glfw.h>
#endif // GLFW2
#endif // INOVESA_USE_OPENGL
#include <vector>

#include "PS/PhaseSpace.hpp"
#include "IO/GUI/GUIElement.hpp"

namespace vfps {

class Display;

class DisplayException : public std::exception {
public:
    DisplayException(std::string msg) : _msg(msg){}

    const char* what() const noexcept override;

private:
    std::string _msg;
};

/**
 * @brief make_display factory function for Display
 * @param gui
 * @param ofname name of result file (will be used for according log file)
 * @param glversion
 * @return pointer to fully initialized Display (may be non-graphical)
 *
 * This factory function will also initialize loging.
 * When no (graphical) display is wanted make_display may
 * just initialize the log, print some status information, or do nothing
 * but to return a nullptr.
 */
std::unique_ptr<Display> make_display(std::string ofname
                                      #ifdef INOVESA_USE_OPENGL
		                      , bool gui
                                      , uint_fast8_t glversion=0
                                      #endif // INOVESA_USE_OPENGL
                                     );

class Display
{
public:
    /**
     * @brief start_time time stamp of the program start
     *
     * start_time will be used in output and log files to display the
     * progress of execution.
     */
    static std::chrono::system_clock::time_point start_time;

    volatile static bool abort;

public:
    Display() = delete;

    Display(const Display&) = delete;
    Display(Display&&) = delete;

    Display& operator=(const Display&) = delete;
    Display& operator=(Display&&) = delete;

    #ifdef INOVESA_USE_OPENGL
    /**
     * @brief Display initializes OpenGL
     * @param glversion
     */
    Display(uint_fast8_t glversion);
    #endif // INOVESA_USE_OPENGL

    /**
     * @brief ~Display() terminats OpenGL (if used)
     */
    ~Display() noexcept;

    #ifdef INOVESA_USE_OPENGL
    void addElement(std::shared_ptr<GUIElement> newitem);
    #endif // INOVESA_USE_OPENGL

    void draw();

    static void printText(std::string txt, float silentTime=0.0f);

    #ifdef INOVESA_USE_OPENGL
    void takeElement(std::shared_ptr<GUIElement> item);
    #endif // INOVESA_USE_OPENGL

    static std::ofstream logfile;


private:
    #ifdef INOVESA_USE_OPENGL
    #if GLFW_VERSION_MAJOR == 3

    GLFWwindow* openWindow(uint_fast8_t glversion);

    /**
     * @brief _window pointer to window struct
     *
     * must not be deleted, is owned by GLFW runtime
     */
    GLFWwindow* _window;
    #endif // GLFW_VERSION_MAJOR == 3

    std::vector<std::shared_ptr<GUIElement>> _item;
    #endif // INOVESA_USE_OPENGL

    /**
     * @brief _lastmessage
     */
    static std::chrono::system_clock::time_point _lastmessage;
};

} // namespace vfps

#endif // DISPLAY_HPP

