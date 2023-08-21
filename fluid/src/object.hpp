#pragma once

#include <vector>
#include <glm/glm.hpp>

class Quad {
    public:
        Quad(glm::vec3 m_position = glm::vec3(0.0, 0.0, 0.0), glm::vec2 m_scale = glm::vec2(1.0, 1.0));
        ~Quad();
        void Draw() const;
        glm::vec3 position() const;
        glm::vec3 scale() const;

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
        glm::vec3 m_position;
        glm::vec2 m_scale;
        glm::vec3 m_orientation;

        unsigned int m_VBO, m_EBO, m_VAO;
};
