#include <iostream>
#include <chrono>
#include "glm/glm.hpp"

#include "Message/Message.hpp"
#include "Physics/Physics.hpp"
#include "Renderer/Renderer.hpp"
#include "Input/Input.hpp"
#include "Console/Console.hpp"

#include "Renderer/utils.hpp"
#include "Utils/utils.hpp"

unsigned int w = 1024;
unsigned int h = 768;
Camera camera(glm::vec3(0.0f,0.0f, 5.0f));
float lastX = w / 2.0f;
float lastY = h / 2.0f;
bool firstMouse = true;
float deltaTime = 0.0f;

#define _GRAPH

int main() {
    MessageBus messageBus;
#ifdef _GRAPH
    Renderer renderer(&messageBus);
    Input input(&messageBus, renderer.window);
#endif
    Console console;
    Fluid fluid("inflow", &messageBus);

    // messageBus.add_receiver(physics.read_message);
#ifdef _GRAPH
    messageBus.add_receiver(renderer.read_message);
    messageBus.add_receiver(input.read_message);
#endif
    messageBus.add_receiver(console.read_message);

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
    quads.push_back(
            {glm::vec3(xr,ys - 0*(size_y + space), -1.0f),
            glm::vec2(2.0f, 1.0f),
            fluid.V}
            );
    quads.push_back(
            {glm::vec3(xl,ys - 1*(size_y + space), -1.0f),
            glm::vec2(2.0f, 1.0f),
            fluid.P}
            );


    bool running = true;
    int n=0;
    while (running) {
        auto startFrameTime = std::chrono::system_clock::now();
#ifdef _GRAPH
        input.update();
        processInput(renderer.window);
#endif

        { // Timer t;
        fluid.update();
        }

#ifdef _GRAPH
        running = renderer.update(quads, fluid.constants.imax+2, fluid.constants.jmax+2);
        glfwPollEvents();
#endif

        messageBus.dispatch();

        auto endFrameTime = std::chrono::system_clock::now();
        auto durationFrame = std::chrono::duration_cast<std::chrono::microseconds>(endFrameTime-startFrameTime).count();
        // std::cout << durationFrame/1000. << " ms\t\r" << std::endl;
        deltaTime = durationFrame/1000000.;
        n++;

    };

    return 0;
}
