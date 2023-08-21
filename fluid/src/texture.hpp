#pragma once
#include <string>

class Texture {
    public:
        Texture();
        Texture(std::string texture_file);
        void use() const;

    private:
        void initialize();
        void read_texture_file(std::string texture);
        int width, height, channels;
        unsigned int textureID;
};
