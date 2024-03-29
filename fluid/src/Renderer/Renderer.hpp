#pragma once
#include "glad/glad.h"
#include <GLFW/glfw3.h>
#include <glm/ext/matrix_clip_space.hpp>
#include <glm/gtx/transform.hpp>
#include <string>
#include <iostream>

#include "Renderer/callbacks.hpp"
#include "Renderer/shader.hpp"
#include "Renderer/texture.hpp"
#include "Renderer/object.hpp"
#include "Renderer/camera.hpp"
#include "Message/Message.hpp"

extern Camera camera;

class Renderer: Sender {
public:
    std::string windowTitle = "Simulation";
    unsigned int w = 1024;
    unsigned int h = 768;

    GLFWwindow* window;

    Renderer(MessageBus *mBus): Sender(mBus) {
        window = initialize_context();
        shader = new Shader(                    \
            "resources/shader/vertex.vs",   \
            "resources/shader/fragment.fs"  \
            );

        texture = new Texture();
        texture->use();

        shader->use();
        shader->setInt("texture1", 0);

        q = new Quad();

        glEnable(GL_DEPTH_TEST);
    }

    static void read_message(Message m){
    };

    ~Renderer() {
        glfwDestroyWindow(window);
        delete texture;
        delete shader;
        delete q;
    };

    bool update(const std::vector<Object>& quads, int width, int height) {
        bool running = !glfwWindowShouldClose(window);
        glm::mat4 proj = glm::perspective(glm::radians(camera.Zoom), (float)w/h, 0.1f, 180.0f);
        glm::mat4 view = camera.GetViewMatrix();
        glm::mat4 model = glm::mat4(1.0f);

        shader->setMat4f("proj", proj);
        shader->setMat4f("view", view);
        shader->setMat4f("model", model);
        glClearColor(0.2f, 0.3f, 0.4f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        for (std::vector<Object>::size_type i=0; i!=quads.size(); ++i) {
            texture->load_texture(quads[i].tensor->data(), width, height, 1);
            model = glm::mat4(1.0f);
            model = glm::scale(model, glm::vec3(quads[i].scale, 1.0f));
            model = glm::translate(model, quads[i].position);
            shader->setMat4f("model", model);
            q->Draw();
        }
        glfwSwapBuffers(window);
        return running;
    };

private:
    Shader *shader;
    Texture *texture;
    Quad *q;

    GLFWwindow* initialize_context() {
        glfwInit();
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef _MAC
        glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

        GLFWwindow *window = glfwCreateWindow(w, h, windowTitle.c_str(), NULL, NULL);

        if (window == nullptr)
        {
            std::cout << "Failed to create GLFW window" << std::endl;
            glfwTerminate();
            return nullptr;
        }
        glfwMakeContextCurrent( window );
        glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
        glfwSetCursorPosCallback(window, mouse_callback);
        glfwSetScrollCallback(window, scroll_callback);

        // glad: load all OpenGL function pointers
        // ---------------------------------------
        if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
        {
            std::cout << "Failed to initialize GLAD" << std::endl;
            return nullptr;
        }
        return window;
    };
};
