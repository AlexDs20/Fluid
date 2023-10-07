#pragma once

struct Parameters {
    // Fluid
    const float Re = 10000;
    const float u0 = 0.0f;
    const float v0 = 0.0f;
    const float p0 = 0.0f;
    const float gx = 0.0;       // not currently used
    const float gy = 0.0;

    // time
    // m
    float t = 0;
    const float t_max = 100;
    float dt = 0.04;
    const float tau = 0.2;

    // grid
    int dims = 2;
    const float length_x = 2.0f;
    const float length_y = 1.0f;
    const int imax = 150;
    const int jmax = 75;
    const float dx = length_x / (float)imax;
    const float dy = length_y / (float)jmax;

    // pressure
    int it = 0;
    int it_max = 10;
    float rit = 0.0f;
    float eps = 0.001f;
    float omega = 1.5;
    float gamma = 0.9f;
    float norm_p0 = 0.5f;
    const float alpha = 0.5f;

    // Other stuff
    int n = 0;
};
