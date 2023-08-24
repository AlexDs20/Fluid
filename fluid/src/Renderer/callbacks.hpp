#pragma once
#include <GLFW/glfw3.h>

// Callback to resize Viewport
void framebuffer_size_callback(GLFWwindow* window, int width, int height);

// Callback for the mouse input
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
