// #include <glad/glad.h>
// #include <GLFW/glfw3.h>

#include <glm/ext/matrix_clip_space.hpp>
#include <iostream>
#include <string>
#include <math.h>

#include "Message/Message.hpp"
#include "Physics/Physics.hpp"
#include "Renderer/Renderer.hpp"

#include "Renderer/setupGL.hpp"
#include "Renderer/utils.hpp"
#include "Renderer/object.hpp"
#include "Renderer/shader.hpp"
#include "Renderer/camera.hpp"
#include "Renderer/texture.hpp"

#include "Math/tensor.hpp"
#include "Physics/calculate.hpp"

#include "parameters.hpp"

#include <chrono>
#include <thread>


unsigned int w = 1920;
unsigned int h = 1080;
Camera camera(glm::vec3(0.0f,0.0f, 5.0f));
float lastX = w / 2.0f;
float lastY = h / 2.0f;
bool firstMouse = true;
float deltaTime = 0.0f;
float lastFrame = 0.0f;


int main() {
    MessageBus messageBus;
    Physics physics(&messageBus);
    Renderer renderer;
    messageBus.add_receiver(physics.read_message);

    Fluid fluid;

    //--------------------
    //  Scene
    std::vector<Placement> quads;
    float xl = -0.55f, xr = 0.55f;
    float ys = 0.55f;
    float size_y = 1.0f;
    float space = 0.25f;
    quads.push_back(Placement({xl,ys - 0*(size_y + space), -1.0f}, {2.0f, 1.0f}));
    quads[0].tensor = fluid.U;

    //--------------------
    //  Renderer
    GLFWwindow* window = renderer.window;

    //--------------------
    //  Scene + entity

    bool running = true;
    while (running) {
        running = !glfwWindowShouldClose(window);
        float currentFrame = glfwGetTime();

        processInput(window);

        fluid.update();

        //------------------------------
        //  RENDERER
        running = renderer.update(quads);
        glfwPollEvents();

        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
    };

    return 0;
}
