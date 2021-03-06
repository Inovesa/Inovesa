// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#pragma once

#include "defines.hpp"

#include <array>
#include <chrono>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <type_traits>

#if INOVESA_ENABLE_INTERRUPT == 1
#include<csignal> // for SIGINT handling
#endif

#if INOVESA_USE_OPENGL == 1

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

#if INOVESA_USE_OPENCL == 1
#include "CL/OpenCLHandler.hpp"
#endif

#include "PS/PhaseSpace.hpp"
#include "IO/GUI/GUIElement.hpp"

namespace vfps {

class Display;

class DisplayException : public std::exception {
public:
    explicit DisplayException(const std::string& msg) : _msg(msg){}

    const char* what() const noexcept override;

private:
    const std::string _msg;
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
std::unique_ptr<Display> make_display( std::string ofname
                                     , bool gui=false
                                     , uint_fast8_t glversion=0
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

    volatile static bool silent_mode;

public:
    Display() = delete;

    Display(const Display&) = delete;
    Display(Display&&) = delete;

    Display& operator=(const Display&) = delete;
    Display& operator=(Display&&) = delete;

    #if INOVESA_USE_OPENGL == 1
    /**
     * @brief Display initializes OpenGL
     * @param glversion
     */
    explicit Display(uint_fast8_t glversion);
    #endif // INOVESA_USE_OPENGL

    /**
     * @brief ~Display() terminats OpenGL (if used)
     */
    ~Display() noexcept;

    #if INOVESA_ENABLE_INTERRUPT == 1
    static void SIGINT_handler(int) {
        Display::abort = true;
    }
    #endif // INOVESA_ENABLE_INTERRUPT

    #if INOVESA_USE_OPENGL == 1
    void addElement(std::shared_ptr<GUIElement> newitem);
    #endif // INOVESA_USE_OPENGL

    void draw();

    /**
     * @brief prints a timestamp followed by a text to both,
     *        stdout and logfile (if present)
     *
     * @param txt text to print
     * @param newline Start new line after printing (true) or print next text on top (false)?
     *        This is just for printing and has no effect for the logfile.
     * @param silentTime skip if last printout was less then silentTime seconds ago.
     *
     * The timestamp is designed to be in seconds after program start.
     *
     * @todo replace by << operator
     */
    static void printText(std::string txt,
                          bool newline = true,
                          float silentTime=0.0f);

    #if INOVESA_USE_OPENGL == 1
    void takeElement(std::shared_ptr<GUIElement> item);
    #endif // INOVESA_USE_OPENGL

    static std::ofstream logfile;


private:
    #if INOVESA_USE_OPENGL == 1
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

