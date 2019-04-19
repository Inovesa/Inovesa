// SPDX-License-Identifier: GPL-3.0-or-later
/*
 * Copyright (c) Patrik Schönfeldt
 * Copyright (c) Karlsruhe Institute of Technology
 */

#include <sstream>

#include "HelperFunctions.hpp"
#include "IO/Display.hpp"


std::unique_ptr<vfps::Display> vfps::make_display(std::string ofname
                                                  #if INOVESA_USE_OPENGL == 1
                                                  , bool gui
                                                  , uint_fast8_t glversion
                                                  #endif // INOVESA_USE_OPENGL
                                                 )
{
    const std::time_t start_ctime
            = std::chrono::system_clock::to_time_t(Display::start_time);
    std::stringstream sstream;
    sstream << std::put_time(std::localtime(&start_ctime),"%FT%T%z");

    std::string timestring = sstream.str();

    if (!ofname.empty()) {
        Display::logfile.open(ofname+".log");
    }
    Display::printText("Started Inovesa ("
                       +vfps::inovesa_version()+") at "+timestring);
    if (!ofname.empty()) {
        Display::printText("Will create log at \""+ofname+".log\".");
    }
    #if INOVESA_USE_OPENGL == 1
    if (gui) {
        return std::make_unique<Display>(glversion);
    }
    #endif // INOVESA_USE_OPENGL
    return nullptr;
}

#if INOVESA_USE_OPENGL == 1
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
    glfwOpenWindow( 512, 512,6,5,6,0,0,0, GLFW_WINDOW);
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
    #if INOVESA_USE_OPENGL == 1
    // Close OpenGL window and terminate GLFW
    glfwTerminate();
    #endif // INOVESA_USE_OPENGL
}

#if INOVESA_USE_OPENGL == 1
void vfps::Display::addElement(std::shared_ptr<GUIElement> newitem)
{
    _item.push_back(newitem);
}
#endif // INOVESA_USE_OPENGL

#if INOVESA_USE_OPENGL == 1
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

void vfps::Display::printText(std::string txt, bool newline, float silentTime)
{
    std::stringstream message;
    std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
    if (std::chrono::duration<float>(now-_lastmessage).count() >= silentTime) {
        message.setf( std::ios::fixed, std:: ios::floatfield );
        message.precision(2);
        message << "[ " << std::setw(9)
                << std::chrono::duration<float>(now-start_time).count()
                << " ]: "
                << txt;
        _lastmessage = now;
        if (newline) {
            std::cout << message.str() << std::endl;
        } else {
            std::cout << message.str() << "\r";
        }
        std::cout.flush();
        if (logfile.is_open()) {
            logfile << message.str() << std::endl;
            logfile.flush();
        }
    }
}

#if INOVESA_USE_OPENGL == 1
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
