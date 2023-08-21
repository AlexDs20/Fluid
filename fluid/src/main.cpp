#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/ext/matrix_clip_space.hpp>
#include <iostream>
#include <string>
#include <math.h>

#include "setupGL.hpp"
#include "utils.hpp"
#include "object.hpp"
#include "shader.hpp"
#include "camera.hpp"
#include "texture.hpp"

#include "tensor.hpp"
#include "functions.hpp"
#include "parameters.hpp"

unsigned int w = 1920;
unsigned int h = 1080;
Camera camera(glm::vec3(0.0f, 0.0f, 3.0f));
float lastX = w / 2.0f;
float lastY = h / 2.0f;
bool firstMouse = true;
float deltaTime = 0.0f;
float lastFrame = 0.0f;

int main(int argc, char** argv) {
    std::string config;
    if (argc>1)
        config = argv[1];

    //--------------------
    std::string window_title("Simulation");
    GLFWwindow* window = setupGL(window_title, w, h);

    Shader shader( \
            "resources/shader/vertex.vs", \
            "resources/shader/fragment.fs" \
    );
    Texture texture(
            "resources/textures/container.jpg" \
    );
    int sx = 32;
    int sy = 32;
    Tensor U({3, sx, sy}, 0.05f);
    for (int i=0; i!=sx;++i)
        for (int j=0; j!=sy;++j) {
            U({0, i, j}) = 0.0;
            U({1, i, j}) = 0.9;
            U({2, i, j}) = 0.0;
        }
    texture.load_texture(U.data(), sx, sy, 1);

    Quad q({0.0f, 0.0f, -1.0f}, {2.0f, 1.0f});
    shader.use();
    shader.setInt("texture1", 0);
    texture.use();

    glm::mat4 proj = glm::perspective(glm::radians(camera.Zoom), (float)w/h, 0.1f, 100.0f);
    glm::mat4 view = camera.GetViewMatrix();
    glm::mat4 model = glm::mat4(1.0f);
    model = glm::scale(model, q.scale());
    model = glm::translate(model, q.position());
    // model = glm::rotate(model, );

    shader.setMat4f("proj", proj);
    shader.setMat4f("view", view);
    shader.setMat4f("model", model);

    glEnable(GL_DEPTH_TEST);

    while (!glfwWindowShouldClose(window)) {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        // std::cout << deltaTime*1000 << " ms\r";
        // std::cout.flush();

        processInput(window);

        glClearColor(0.2f, 0.3f, 0.4f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        q.Draw();

        glfwSwapBuffers(window);
        glfwPollEvents();
    };
    glfwTerminate();
    //--------------------

    // Parameters p;

    // Initialize the fields
    // TODO: Investigate if order of the dim / grid or if grouping the fields together changes the perf.
    // Tensor U({p.imax+2, p.jmax+2}, p.u0);
    // Tensor V({p.imax+2, p.jmax+2}, p.v0);
    // Tensor P({p.imax+2, p.jmax+2}, p.p0);
    // Tensor EXT_X({p.imax+2, p.jmax+2});     // External forces
    // Tensor EXT_Y({p.imax+2, p.jmax+2});     // External forces
    // Tensor F({p.imax+2, p.jmax+2});
    // Tensor G({p.imax+2, p.jmax+2});
    // Tensor RHS({p.imax+2, p.jmax+2});
    // Tensor OBS({p.imax+2, p.jmax+2});

    // while (p.t < p.t_max) {
    //     std::cout << p.t << "/" << p.t_max << "\t" << p.dt << std::endl;
    //     p.dt = adaptive_time_step_size(U, V, p.dx, p.dy, p.Re, p.tau, p.dt);
    //     set_boundary_values(U, V);
    //     set_object_boundary_values(U, V, OBS);
    //     compute_FG(F, G, U, V, p.dt, p.Re, p.dx, p.dy, p.gamma);
    //     compute_rhs_pressure(RHS, F, G, p.dx, p.dy, p.dt);

    //     p.it = 0;
    //     while (p.it<p.it_max && p.rit > p.eps * p.norm_p0) {
    //         SOR(P, p.rit, RHS, p.omega, p.dx, p.dy);
    //         ++p.it;
    //     }

    //     compute_uv(U, V, F, G, P, p.dx, p.dy, p.dt);

    //     p.t += p.dt;
    //     p.n += 1;
    // }

    //--------------------
    glfwDestroyWindow(window);
    //--------------------

    return 0;
}
