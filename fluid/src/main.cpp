#include <iostream>
#include <glm/glm.hpp>

#include "Message/Message.hpp"
#include "Physics/Physics.hpp"
#include "Renderer/Renderer.hpp"
#include "Input/Input.hpp"
#include "Console/Console.hpp"

#include "Renderer/utils.hpp"

#include <chrono>


unsigned int w = 1920;
unsigned int h = 1080;
Camera camera(glm::vec3(0.0f,0.0f, 5.0f));
float lastX = w / 2.0f;
float lastY = h / 2.0f;
bool firstMouse = true;
float deltaTime = 0.0f;


int main() {
    MessageBus messageBus;
    Physics physics(&messageBus);
    // Renderer renderer(&messageBus);
    // Input input(&messageBus, renderer.window);
    Console console;

    messageBus.add_receiver(physics.read_message);
    // messageBus.add_receiver(renderer.read_message);
    // messageBus.add_receiver(input.read_message);
    messageBus.add_receiver(console.read_message);

    Fluid fluid;

    //--------------------
    //  Scene
    std::vector<Object> quads;
    float xl = -0.55f, xr = 0.55f;
    float ys = 0.55f;
    float size_y = 1.0f;
    float space = 0.25f;

    quads.push_back(
            {glm::vec3(xl,ys - 0*(size_y + space), -1.0f),
            glm::vec2(2.0f, 1.0f),
            fluid.U}
            );


    bool running = true;
    while (running) {
        auto startFrameTime = std::chrono::system_clock::now();
        // input.update();
        // processInput(renderer.window);

        fluid.update();

        // running = renderer.update(quads);
        glfwPollEvents();

        messageBus.dispatch();

        auto endFrameTime = std::chrono::system_clock::now();
        auto durationFrame = std::chrono::duration_cast<std::chrono::microseconds>(endFrameTime-startFrameTime).count();
        // std::cout << durationFrame/1000. << " ms\t\r" << std::endl;
        deltaTime = durationFrame/1000000.;
    };

    return 0;
}
