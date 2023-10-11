#include <glad/glad.h>
#include <GLFW/glfw3.h>

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
Camera camera(glm::vec3(0.0f,-1.5f, 7.0f));
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


    //--------------------
    //  RENDERER
    //--------------------
    GLFWwindow* window = renderer.window;

    Shader shader( \
            "resources/shader/vertex.vs", \
            "resources/shader/fragment.fs" \
    );

    Texture texture;

    std::vector<Placement> quads;
    float xl = -0.55f, xr = 0.55f;
    float ys = 0.55f;
    float size_y = 1.0f;
    float space = 0.25f;

    quads.push_back(Placement({xl,ys - 0*(size_y + space), -1.0f}, {2.0f, 1.0f}));
    quads.push_back(Placement({xr,ys - 0*(size_y + space), -1.0f}, {2.0f, 1.0f}));
    quads.push_back(Placement({xl,ys - 1*(size_y + space), -1.0f}, {2.0f, 1.0f}));

    Quad q;

    shader.use();
    shader.setInt("texture1", 0);
    texture.use();

    glEnable(GL_DEPTH_TEST);


    Fluid fluid;
    Parameters p;
    quads[0].tensor = fluid.U;
    quads[1].tensor = fluid.V;
    quads[2].tensor = fluid.P;

    float dt;
    int n = 0;
    float simulated_t = 0;
    float elapsed_t = 0;
    while (!glfwWindowShouldClose(window) & (simulated_t<p.t_max)) {
        float currentFrame = glfwGetTime();
        std::cout << "TIME: " << "Elapsed " << elapsed_t << "\t" << "Simulated: " << simulated_t << "\t";
        std::cout << "frame: " << deltaTime*1000 << " ms\t" << "simulated step: " << dt << "\t" << "iter: " << n << "\t" << std::endl;
        std::cout.flush();

        //------------------------------
        //  RENDERER
        processInput(window);
        glm::mat4 proj = glm::perspective(glm::radians(camera.Zoom), (float)w/h, 0.1f, 180.0f);
        glm::mat4 view = camera.GetViewMatrix();
        glm::mat4 model = glm::mat4(1.0f);

        shader.setMat4f("proj", proj);
        shader.setMat4f("view", view);
        shader.setMat4f("model", model);

        fluid.update();

        //------------------------------
        //  RENDERER
        glClearColor(0.2f, 0.3f, 0.4f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        for (std::vector<Placement>::size_type i=0; i!=quads.size(); ++i) {
            texture.load_texture(quads[i].tensor->data(), p.imax+2, p.jmax+2, 1);
            model = glm::mat4(1.0f);
            model = glm::scale(model, glm::vec3(quads[i].scale, 1.0f));
            model = glm::translate(model, quads[i].position);
            shader.setMat4f("model", model);
            q.Draw();
        }
        glfwSwapBuffers(window);
        glfwPollEvents();

        if (n == 1)
        {
            std::chrono::seconds duration(1);
            std::this_thread::sleep_for(duration);
        }
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        elapsed_t += deltaTime;

    };
    glfwDestroyWindow(window);

    return 0;
}
