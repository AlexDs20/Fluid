#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#include <glad/glad.h>
#include <string>
#include <iostream>
#include "texture.hpp"

Texture::Texture(std::string texture) {
    initialize();
    read_texture_file(texture);
};

Texture::Texture() {
    initialize();
};

void Texture::load_texture(const std::vector<float>& data, int width, int height, int channels) {
    this->width = width;
    this->height = height;
    this->channels = channels;

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, width, height, 0, GL_RGB, GL_FLOAT, data.data());
    glGenerateMipmap(GL_TEXTURE_2D);
}

void Texture::use() const {
    // TODO: This TEXTURE0 should be more flexible
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, textureID);
}

void Texture::read_texture_file(std::string texture) {
    unsigned char* data = stbi_load(texture.c_str(), &width, &height, &channels, 0);

    if (data) {
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);
    } else {
        std::cout << "Failed to load texture!" << std::endl;
    }
    stbi_image_free(data);
};

void Texture::initialize() {
    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_2D, textureID);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
};
