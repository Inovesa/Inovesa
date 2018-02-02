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

#include <sstream>

#include "MessageStrings.hpp"
#include "IO/Display.hpp"


std::unique_ptr<vfps::Display> vfps::make_display(std::string ofname
                                                  #ifdef INOVESA_USE_OPENGL
                                                  , bool gui
                                                  , uint_fast8_t glversion
                                                  #endif // INOVESA_USE_OPENGL
						  )
{
    const std::time_t start_ctime
            = std::chrono::system_clock::to_time_t(Display::start_time);
    std::stringstream sstream;
    sstream << std::ctime(&start_ctime);

    std::string timestring = sstream.str();
    timestring.resize(timestring.size()-1);

    if (!ofname.empty()) {
        Display::logfile.open(ofname+".log");
    }
    Display::printText("Started Inovesa ("
                       +vfps::inovesa_version()+") at "+timestring);
    if (!ofname.empty()) {
        Display::printText("Will create log at \""+ofname+".log\".");
    }
    #ifdef INOVESA_USE_OPENGL
    if (gui) {
        return std::make_unique<Display>(glversion);
    }
    #endif // INOVESA_USE_OPENGL
    return nullptr;
}

#ifdef INOVESA_USE_OPENGL
vfps::Display::Display(uint_fast8_t glversion)
    #if GLFW_VERSION_MAJOR == 3
        :
        _window(nullptr)
    #endif // GLFW_VERSION_MAJOR == 3
{
    if (glversion == 0) {
        #if defined(__APPLE__) || defined(__MACOSX)
        glversion = 2;
        #else
        glversion = 3;
        #endif
    }
    glfwInit();

    // Open a window and create its OpenGL context
    #if GLFW_VERSION_MAJOR < 3
    glfwOpenWindowHint(GLFW_FSAA_SAMPLES, 4);
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MAJOR, 2);
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 1);
    glfwOpenWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwOpenWindow( 512, 512,6,5,6,0,0,0, GLFW_WINDOW);GLFWwindow
    glfwSetWindowTitle("Inovesa");
    GUIElement::glversion = 2;
    #else // GLFW3
    _window = openWindow(glversion);
    if( _window == nullptr ) {
        glfwTerminate();
        throw DisplayException("Failed to initialize GLFW.");
        return;
    }
    glfwMakeContextCurrent(_window);
    #endif // GLFW3

    // Initialize GLEW
    glewExperimental = true; // Needed for core profile
    if (glewInit() != GLEW_OK) {
        throw DisplayException("Failed to initialize GLEW");
    }

    // Ensure we can capture the escape key being pressed below
    #if GLFW_VERSION_MAJOR < 3
    glfwEnable(GLFW_STICKY_KEYS);
    #else
    glfwSetInputMode(_window, GLFW_STICKY_KEYS, GL_TRUE);
    #endif // GLFW_VERSION_MAJOR == 3

    // White background
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f);

    // Enable depth test
    glEnable(GL_DEPTH_TEST);
    // Accept fragment if it closer to the camera than the former one
    glDepthFunc(GL_LESS);
}
#endif // INOVESA_USE_OPENGL

vfps::Display::~Display() noexcept
{
    #ifdef INOVESA_USE_OPENGL
    // Close OpenGL window and terminate GLFW
    glfwTerminate();
    #endif // INOVESA_USE_OPENGL
}

#ifdef INOVESA_USE_OPENGL
void vfps::Display::addElement(std::shared_ptr<GUIElement> newitem)
{
    _item.push_back(newitem);
}
#endif // INOVESA_USE_OPENGL

#ifdef INOVESA_USE_OPENGL
void vfps::Display::draw() {
    if (! glfwWindowShouldClose(_window)) {
        // Clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // draw GUIElements
        for (auto i : _item) {
            i->draw();
        }

        // Swap buffers
        #if GLFW_VERSION_MAJOR == 3
        glfwSwapBuffers(_window);
        #else
        glfwSwapBuffers();
        #endif
        glfwPollEvents();
    } else {
        Display::abort = true;
    }
}
#endif // INOVESA_USE_OPENGL

void vfps::Display::printText(std::string txt, float silentTime)
{
    std::stringstream message;
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    if (std::chrono::duration<float>(now-_lastmessage).count() > silentTime) {
        message.setf( std::ios::fixed, std:: ios::floatfield );
        message.precision(2);
        message << "[ " << std::setw(9)
                << std::chrono::duration<float>(now-start_time).count()
                << " ]: "
                << txt
                << std::endl;
        _lastmessage = now;
    }
    std::cout << message.str();
    std::cout.flush();
    if (logfile.is_open()) {
        logfile << message.str();
        logfile.flush();
    }
}

#ifdef INOVESA_USE_OPENGL
void vfps::Display::takeElement(std::shared_ptr<GUIElement> item)
{
    for (size_t i=0; i< _item.size(); i++) {
        if (_item[i] == item) {
            _item.erase(_item.begin()+i);
            i--;
        }
    }
}

GLFWwindow* vfps::Display::openWindow(uint_fast8_t glversion)
{
    std::string title("Inovesa");
    GUIElement::glversion = glversion;
    switch (glversion) {
    case 2:
        glfwWindowHint(GLFW_SAMPLES, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
        break;
    case 3:
    default:
        glfwWindowHint(GLFW_SAMPLES, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
        title+=" (GL3)";
        break;
    }
    return glfwCreateWindow( 512, 512, title.c_str(), NULL, NULL);
}
#endif // INOVESA_USE_OPENGL

std::chrono::system_clock::time_point vfps::Display::start_time;

volatile bool vfps::Display::abort(false);

std::chrono::system_clock::time_point vfps::Display::_lastmessage;

std::ofstream vfps::Display::logfile;

const char*
vfps::DisplayException::what() const noexcept
    { return _msg.c_str(); }
