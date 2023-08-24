#pragma once
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include "Renderer/camera.hpp"

extern Camera camera;
extern float lastX, lastY;
extern bool firstMouse;

extern float deltaTime;
extern float lastFrame;

// Process Input and how to handle it
void processInput(GLFWwindow* window);

