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


int main(int argc, char** argv) {
    std::string config;
    if (argc>1)
        config = argv[1];

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
    // quads.push_back(Placement({xr,ys - 1*(size_y + space), -1.0f}, {2.0f, 1.0f}));
    // quads.push_back(Placement({xl,ys - 2*(size_y + space), -1.0f}, {2.0f, 1.0f}));
    // quads.push_back(Placement({xr,ys - 2*(size_y + space), -1.0f}, {2.0f, 1.0f}));
    // quads.push_back(Placement({xl,ys - 3*(size_y + space), -1.0f}, {2.0f, 1.0f}));
    // quads.push_back(Placement({xr,ys - 3*(size_y + space), -1.0f}, {2.0f, 1.0f}));
    // quads.push_back(Placement({xl,ys - 4*(size_y + space), -1.0f}, {2.0f, 1.0f}));
    // quads.push_back(Placement({xr,ys - 4*(size_y + space), -1.0f}, {2.0f, 1.0f}));
    // quads.push_back(Placement({xl,ys - 5*(size_y + space), -1.0f}, {2.0f, 1.0f}));

    Quad q;

    shader.use();
    shader.setInt("texture1", 0);
    texture.use();

    glEnable(GL_DEPTH_TEST);

    //--------------------
    // Physics physics;
    Parameters p;
    // Initialize the fields
    Tensor U({p.imax+2, p.jmax+2}, p.u0);
    Tensor V({p.imax+2, p.jmax+2}, p.v0);
    Tensor P({p.imax+2, p.jmax+2}, p.p0);
    Tensor F({p.imax+2, p.jmax+2}, p.u0);
    Tensor G({p.imax+2, p.jmax+2}, p.v0);
    Tensor RHS({p.imax+2, p.jmax+2}, 0.0f);

    // -----------
    // TODO
    Boundary boundary = {
        Tensor({p.imax+2, p.jmax+2}, 0.0f),         // cellType: 0 -> fluid
        Tensor({p.imax+2, p.jmax+2}, 0.0f),
        Tensor({p.imax+2, p.jmax+2}, 0.0f),
        Tensor({p.imax+2, p.jmax+2}, 0.0f),
        Tensor({p.imax+2, p.jmax+2}, 0.0f),
        Tensor({p.imax+2, p.jmax+2}, 0.0f),
        Tensor({p.imax+2, p.jmax+2}, 0.0f),
    };

    quads[0].tensor = &U;
    quads[1].tensor = &V;
    quads[2].tensor = &P;
    // quads[3].tensor = &boundary.cellType;
    // quads[4].tensor = &F;
    // quads[5].tensor = &G;
    // quads[6].tensor = &RHS;
    // quads[7].tensor = &boundary.NCell;
    // quads[8].tensor = &boundary.SCell;
    // quads[9].tensor = &boundary.WCell;
    // quads[10].tensor = &boundary.ECell;

    float dt;

    // Set flags for edges and obstacle
    set_constant_flags(boundary);

    while (!glfwWindowShouldClose(window) & (p.t<p.t_max)) {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        std::cout << deltaTime*1000 << " ms\t" << p.t << " " << "\t" << dt << "\t" << p.n << std::endl;
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
        // physics.update(dt);
        // set_boundary_flags(boundary);
        set_boundary_values(U, V, boundary);
        set_specific_boundary_values(U, V);

        dt = adaptive_time_step_size(U, V, p.dx, p.dy, p.Re, p.tau, p.dt);

        compute_FG(F, G, U, V, dt, p.Re, p.dx, p.dy, p.gamma, boundary);

        compute_rhs_pressure(RHS, F, G, p.dx, p.dy, dt, boundary);

        int it = 0;
        float rit = 0;
        do {
            ++it;
            SOR(P, rit, RHS, p.omega, p.dx, p.dy, boundary);
        } while (it<5);
        // } while (it<p.it_max && rit > p.eps * p.norm_p0);

        compute_uv(U, V, F, G, P, p.dx, p.dy, dt, boundary);

        p.t += p.dt;
        p.n++;
        //------------------------------


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

        if (p.n == 1)
        {
            std::chrono::seconds duration(1);
            std::this_thread::sleep_for(duration);
        }

    };
    // glfwTerminate();
    glfwDestroyWindow(window);

    return 0;
}

