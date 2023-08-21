#pragma once

struct Parameters {
    // Fluid
    float Re = 100.0f;
    float u0 = 0.0f;
    float v0 = 0.0f;
    float p0 = 0.0f;

    // time
    float t = 0;
    float t_max = 100;
    float dt=0.5;
    float tau = 0.5;

    // grid
    int dims = 2;
    int length_x = 2;
    int length_y = 1;
    int imax = 10;
    int jmax = 5;
    float dx = (float)length_x / imax;
    float dy = (float)length_y / jmax;

    // pressure
    int it = 0;
    int it_max = 10;
    float rit = 0.0f;
    float eps = 0.01f;
    float omega = 0.7f;
    float gamma = 0.5f;
    float norm_p0 = 0.5f;

    // Other stuff
    int n = 0;
};
