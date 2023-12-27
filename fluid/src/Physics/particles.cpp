#include <math.h>
#include <vector>
#include "Math/Matrix.hpp"


inline float
bilinear_interp(float x, float y, float dxinv, float dyinv,                           \
        float xMin, float xMax, float yMin, float yMax, \
        float fieldLowLeft, float fieldLowRight, float fieldUpLeft, float fieldUpRight) {

    float out = dxinv*dyinv * (                                       \
          (xMax - x) * (yMax - y) * fieldLowLeft + (x - xMin) * (yMax - y) * fieldLowRight + \
          (xMax - x) * (y - yMin) * fieldUpLeft  + (x - xMin) * (y - yMin) * fieldUpRight);

    return out;
}

inline float
interpolate_field_U(float& x, float& y, const float dx, const float dy, const Matrix& U) {
    float out;
    static const float dxinv = 1 / dx;
    static const float dyinv = 1 / dy;

    const int i = (int)(x * dxinv) + 1;
    const int j = (int)( (y+dy*0.5f) * dyinv) + 1;

    const float x1 = (i-1)*dx;
    const float x2 = i*dx;
    const float y1 = (j-3.5f)*dy;
    const float y2 = (j-0.5f)*dy;

    const float u1 = U(i-1, j-1);
    const float u2 = U(i, j-1);
    const float u3 = U(i-1, j);
    const float u4 = U(i, j);

    out = bilinear_interp(x, y, dxinv, dyinv, x1, x2, y1, y2, u1, u2, u3, u4);

    return out;
}


inline float
interpolate_field_V(float& x, float& y, const float dx, const float dy, const Matrix& V) {
    float out;
    static const float dxinv = 1 / dx;
    static const float dyinv = 1 / dy;

    const int i = (int)( (x + dx*0.5f) * dxinv) + 1;
    const int j = (int)(y * dyinv) + 1;

    const float x1 = (i-3.5f)*dx;
    const float x2 = (i-0.5f)*dx;
    const float y1 = (j-1)*dy;
    const float y2 = j*dy;

    const float v1 = V(i-1, j-1);
    const float v2 = V(i, j-1);
    const float v3 = V(i-1, j);
    const float v4 = V(i, j);

    out = bilinear_interp(x, y, dxinv, dyinv, x1, x2, y1, y2, v1, v2, v3, v4);

    return out;
}


inline void
forward_Euler(float& x, const float dt, const float u) {
    x += dt * u;
}


inline void
forward_Euler(float& x, float& y, const float dt, const float dx, const float dy, const Matrix& U, const Matrix& V) {
    float u = interpolate_field_U(x, y, dx, dy, U);
    float v = interpolate_field_V(x, y, dx, dy, U);

    forward_Euler(x, dt, u);
    forward_Euler(y, dt, v);
}


// This will contain all the particles current position
// x[i], y[i] is position of particle i
struct Particles {
    float xs, ys, xe, ye;
    float xDir, yDir;
    std::vector<float> x;
    std::vector<float> y;
};


inline float
set_pos(const int i, const float start, const float dir) {
    return start + (i+1) * dir;
}


// Initialize N particles linearly distributed between Start and End position, not included
Particles
init_particles(const int N, const float xStart, const float yStart, const float xEnd, const float yEnd) {
    const float segmentLength = sqrt( pow(yEnd-yStart, 2) + pow(xEnd-xStart, 2) );
    const float incrementLength = segmentLength / (N+1);
    const float xDir = incrementLength * (xEnd - xStart) / segmentLength;
    const float yDir = incrementLength * (yEnd - yStart) / segmentLength;

    Particles particles{xStart, yStart, xEnd, yEnd, xDir, yDir, {}, {}};
    particles.x.reserve(N);
    particles.y.reserve(N);

    for (int i=0; i<N; i++) {
        particles.x[i] = set_pos(i, xStart, xDir);
        particles.y[i] = set_pos(i, yStart, yDir);
    }
    return particles;
}


void
move_particles(Particles& particles,        \
        const int imax, const int jmax, const int dx, const int dy, const int dt,   \
        const Matrix* U, const Matrix* V) {

    const int imp1 = imax+1;
    const int jmp1 = jmax+1;

    for (std::vector<float>::size_type i=0; i<particles.x.size(); i++) {

        forward_Euler(particles.x[i], particles.y[i], dt, dx, dy, *U, *V);

        if ( ((particles.x[i] > imp1*dx) || (particles.x[i]<0)) ||
             ((particles.y[i] > jmp1*dy) || (particles.y[i]<0)) ) {
            particles.x[i] = set_pos(i, particles.xs, particles.xDir);
            particles.y[i] = set_pos(i, particles.ys, particles.yDir);
        }
    }
}

void particle_tracing();
void streaklines();
