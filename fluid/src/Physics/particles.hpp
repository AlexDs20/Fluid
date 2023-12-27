#pragma once

#include <math.h>
#include <vector>
#include "Math/Matrix.hpp"


inline float
bilinear_interp(float x, float y, float dxinv, float dyinv,                           \
        float xMin, float xMax, float yMin, float yMax, \
        float fieldLowLeft, float fieldLowRight, float fieldUpLeft, float fieldUpRight);

inline float
interpolate_field_U(float& x, float& y, const float dx, const float dy, const Matrix& U);


inline float
interpolate_field_V(float& x, float& y, const float dx, const float dy, const Matrix& V);


inline void
forward_Euler(float& x, const float dt, const float u);


inline void
forward_Euler(float& x, float& y, const float dt, const float dx, const float dy, const Matrix& U, const Matrix& V);


// This will contain all the particles current position
// x[i], y[i] is position of particle i
struct Particles {
    float xs, ys, xe, ye;
    float xDir, yDir;
    std::vector<float> x;
    std::vector<float> y;
};


inline float
set_pos(const int i, const float start, const float dir);


// Initialize N particles linearly distributed between Start and End position, not included
Particles
init_particles(const int N, const float xStart, const float yStart, const float xEnd, const float yEnd);


void
move_particles(Particles& particles,        \
        const int imax, const int jmax, const int dx, const int dy, const int dt,   \
        const Matrix* U, const Matrix* V);

void particle_tracing();
void streaklines();
