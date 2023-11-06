#pragma once
#include <vector>
#include <glm/glm.hpp>
#include "Math/Matrix.hpp"

struct Object {
    glm::vec3 position;
    glm::vec2 scale;
    Matrix* tensor = nullptr;
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
