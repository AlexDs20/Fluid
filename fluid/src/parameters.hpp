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
    const float t_max = 100;
    const float dt_max = 0.04;
    const float tau = 0.3;

    // grid
    const float length_x = 2.0f;
    const float length_y = 1.0f;
    const int imax = 150;
    const int jmax = 75;
    const float dx = length_x / (float)imax;
    const float dy = length_y / (float)jmax;

    // pressure
    const int it_max = 5;
    const float eps = 0.01f;
    const float omega = 1.5;

    // Derivative scheme
    const float gamma = 0.9f;
};
