#pragma once
#include <string>
#include "glad/glad.h"

class Texture {
    public:
        Texture();
        ~Texture();
        Texture(std::string texture_file);
        void load_texture(const float* data, int width, int height, int channels);
        void use() const;

    private:
        void initialize();
        void read_texture_file(std::string texture);
        int width, height, channels;
        GLuint textureID;
        // float *D = nullptr;
};
