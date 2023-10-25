#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include <string>
#include <iostream>

#include "Renderer/texture.hpp"

Texture::Texture(std::string texture) {
    initialize();
    read_texture_file(texture);
};

Texture::Texture() {
    initialize();
};

Texture::~Texture() {
    // delete[] D;
};

void Texture::load_texture(const float* data, int width, int height, int channels) {
    this->width = width;
    this->height = height;
    this->channels = channels;

    // if (D == nullptr) {
    //     D = new float[width*height];
    // }

    // int block_size=64;
    // int b2 = block_size * block_size;
    // int kimax = width / block_size;
    // for (int j=0; j<height; j++) {
    //     for (int i=0; i<width; i++){
    //         int ki = i / block_size;
    //         int kj = j / block_size;
    //         int bi = i - ki * block_size;
    //         int bj = j - kj * block_size;
    //         int index = (ki + kj*kimax)*b2 + bi + bj*block_size;
    //         D[i + j*width] = data[index];
    //     }
    // }

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, width, height, 0, GL_RED, GL_FLOAT, data);
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
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
};
