#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/ext/matrix_clip_space.hpp>
#include <iostream>
#include <string>
#include <math.h>

#include "Renderer/setupGL.hpp"
#include "Renderer/utils.hpp"
#include "Renderer/object.hpp"
#include "Renderer/shader.hpp"
#include "Renderer/camera.hpp"
#include "Renderer/texture.hpp"
#include "Renderer/Renderer.hpp"

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
    //--------------------
    //  RENDERER
    //--------------------
    Renderer renderer;
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

    //--------------------
    //  PHYSICS
    //--------------------
    Parameters p;
    // Initialize the fields
    Tensor U({p.imax+2, p.jmax+2}, p.u0);
    Tensor V({p.imax+2, p.jmax+2}, p.v0);
    Tensor P({p.imax+2, p.jmax+2}, p.p0);
    Tensor F({p.imax+2, p.jmax+2}, 0.0f);
    Tensor G({p.imax+2, p.jmax+2}, 0.0f);
    Tensor RHS({p.imax+2, p.jmax+2}, 0.0f);
    Domain domain({p.imax+2, p.jmax+2}, Cell());

    // Set flags for edges and obstacle     (because obstacles don't move)
    set_constant_flags(domain, p.imax, p.jmax);

    quads[0].tensor = &U;
    quads[1].tensor = &V;
    quads[2].tensor = &P;

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

        //------------------------------
        //  PHYSICS
        set_boundary_values(U, V, domain, p.imax, p.jmax);
        set_specific_boundary_values(U, V, p.imax, p.jmax);
        dt = adaptive_time_step_size(U, V, p.dx, p.dy, p.Re, p.tau, p.dt_max, p.imax, p.jmax);
        compute_FG(F, G, U, V, domain, dt, p.Re, p.dx, p.dy, p.gamma, p.imax, p.jmax, p.gx, p.gy);
        compute_rhs_pressure(RHS, F, G, domain, p.dx, p.dy, dt, p.imax, p.jmax);

        int it = 0;
        float rit = 0;
        do {
            ++it;
            SOR(P, RHS, domain, rit, p.omega, p.dx, p.dy, p.imax, p.jmax);
        } while (it < p.it_max && rit > p.eps);

        compute_uv(U, V, F, G, P, domain, p.dx, p.dy, dt, p.imax, p.jmax);

        simulated_t += dt;
        n++;
        //------------------------------

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

