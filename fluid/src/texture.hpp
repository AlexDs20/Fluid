#pragma once
#include <string>
#include <vector>

class Texture {
    public:
        Texture();
        Texture(std::string texture_file);
        void load_texture(const std::vector<float>& data, int width, int height, int channels);
        void use() const;

    private:
        void initialize();
        void read_texture_file(std::string texture);
        int width, height, channels;
        unsigned int textureID;
};
