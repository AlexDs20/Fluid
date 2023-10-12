#pragma once
#include <vector>
#include <glm/glm.hpp>
#include "Math/tensor.hpp"

struct Placement {
    Placement(glm::vec3 pos = glm::vec3(0.0f, 0.0f, 0.0f), \
            glm::vec2 scale = glm::vec2(1.0f, 1.0f), \
            glm::vec2 orientation = glm::vec2(0.0f, 0.0f)): \
        position(pos), scale(scale), orientation(orientation) {};
    glm::vec3 position;
    glm::vec2 scale;
    glm::vec2 orientation;
    Tensor* tensor;
};

class Quad {
    public:
        Quad();
        ~Quad();
        void Draw() const;

    private:
        void m_initialize();
        const std::vector<float> m_vertices = {
        // POSX     POSY   Texture X Y
            -0.5,   0.5,    0.0, 1.0,
             0.5,   0.5,    1.0, 1.0,
             0.5,  -0.5,    1.0, 0.0,
            -0.5,  -0.5,    0.0, 0.0,
        };
        const std::vector<unsigned int> m_indices = {
            0, 1, 2,
            2, 3, 0,
        };

        unsigned int m_VBO, m_EBO, m_VAO;
};
