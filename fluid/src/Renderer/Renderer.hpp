#pragma once
#include <glad/glad.h>
#include <GLFW//glfw3.h>
#include <string>
#include <iostream>

#include "Renderer/callbacks.hpp"

class Renderer {
    public:
        Renderer() {
        }
    private:

        std::string windowTitle = "Simulation";
        unsigned int w = 1920;
        unsigned int h = 1080;

        GLFWwindow* initialize_context() {
            // ------------------------------------
            //      Init glfw, glad, window
            // ------------------------------------
            // glfw: initialize and configure
            // ----------------------------------
            glfwInit();
            glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
            glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
            glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef _MAC
            glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

            // glfw window creation
            // --------------------
            GLFWwindow *window = glfwCreateWindow(w, h, windowTitle.c_str(), NULL, NULL);

            if (window == NULL)
            {
                std::cout << "Failed to create GLFW window" << std::endl;
                glfwTerminate();
                return nullptr;
            }
            glfwMakeContextCurrent( window );
            glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
            glfwSetCursorPosCallback(window, mouse_callback);
            glfwSetScrollCallback(window, scroll_callback);

            // glad: load all OpenGL function pointers
            // ---------------------------------------
            if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
            {
                std::cout << "Failed to initialize GLAD" << std::endl;
                return nullptr;
            }
            return window;
        };
};
