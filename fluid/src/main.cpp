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
#include "calculate.hpp"
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

    Quad qu({-0.55f, 0.55f, -1.0f}, {2.0f, 1.0f});
    Quad qv({ 0.55f, 0.55f, -1.0f}, {2.0f, 1.0f});
    Quad qp({ 0.00f,-0.55f, -1.0f}, {2.0f, 1.0f});
    shader.use();
    shader.setInt("texture1", 0);
    texture.use();

    glEnable(GL_DEPTH_TEST);

    //--------------------
    Parameters p;
    // Initialize the fields
    Tensor U({p.imax+2, p.jmax+2}, p.u0);
    Tensor V({p.imax+2, p.jmax+2}, p.v0);
    Tensor P({p.imax+2, p.jmax+2}, p.p0);
    Tensor EXT_X({p.imax+2, p.jmax+2});     // External forces
    Tensor EXT_Y({p.imax+2, p.jmax+2});     // External forces
    Tensor F({p.imax+2, p.jmax+2});
    Tensor G({p.imax+2, p.jmax+2});
    Tensor RHS({p.imax+2, p.jmax+2});
    Tensor OBS({p.imax+2, p.jmax+2});


    while (!glfwWindowShouldClose(window) & (p.t<p.t_max)) {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        std::cout << deltaTime*1000 << " ms\t" << p.t << " " << std::endl;
        std::cout.flush();

        processInput(window);
        glm::mat4 proj = glm::perspective(glm::radians(camera.Zoom), (float)w/h, 0.1f, 180.0f);
        glm::mat4 view = camera.GetViewMatrix();
        glm::mat4 model = glm::mat4(1.0f);
        // model = glm::rotate(model, );

        shader.setMat4f("proj", proj);
        shader.setMat4f("view", view);
        shader.setMat4f("model", model);

        //------------------------------
        p.dt = adaptive_time_step_size(U, V, p.dx, p.dy, p.Re, p.tau, p.dt);
        set_boundary_values(U, V);
        // set_object_boundary_values(U, V, OBS);
        compute_FG(F, G, U, V, p.dt, p.Re, p.dx, p.dy, p.gamma);
        compute_rhs_pressure(RHS, F, G, p.dx, p.dy, p.dt);

        p.it = 0;
        while (p.it<p.it_max && p.rit > p.eps * p.norm_p0) {
            SOR(P, p.rit, RHS, p.omega, p.dx, p.dy);
            ++p.it;
        }

        compute_uv(U, V, F, G, P, p.dx, p.dy, p.dt);

        p.t += p.dt;
        p.n += 1;
        //------------------------------


        glClearColor(0.2f, 0.3f, 0.4f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        texture.load_texture(U.data(), p.imax+2, p.jmax+2, 1);
        model = glm::mat4(1.0f);
        model = glm::scale(model, qu.scale());
        model = glm::translate(model, qu.position());
        shader.setMat4f("model", model);
        qu.Draw();

        texture.load_texture(V.data(), p.imax+2, p.jmax+2, 1);
        model = glm::mat4(1.0f);
        model = glm::scale(model, qv.scale());
        model = glm::translate(model, qv.position());
        shader.setMat4f("model", model);
        qv.Draw();

        texture.load_texture(P.data(), p.imax+2, p.jmax+2, 1);
        model = glm::mat4(1.0f);
        model = glm::scale(model, qp.scale());
        model = glm::translate(model, qp.position());
        shader.setMat4f("model", model);
        qp.Draw();

        glfwSwapBuffers(window);
        glfwPollEvents();
    };
    // glfwTerminate();
    glfwDestroyWindow(window);

    return 0;
}
