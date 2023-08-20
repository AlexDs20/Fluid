#include <iostream>
#include <math.h>
#include "tensor.hpp"

float adaptive_time_step_size(const Tensor& U, const Tensor& V, float dx, float dy, float Re, float tau, float dt=0.5f) {
    float ret = std::min(dt,
            (Re / 2.0f) * 1.0f / (1.0f/(dx*dx) + 1.0f/(dy*dy))
            );

    float dt_v_max = std::min(dx / U.amax(), dy / V.amax());

    return tau * std::min(ret, dt_v_max);
};

void set_boundary_values(Tensor& U, Tensor& V, const Tensor& OBS){
    for (int i=0; i!=OBS.shape(0); ++i) {
        for (int j=0; j!=OBS.shape(1); ++j) {
        }
    }
};

void compute_FG(){};

void compute_rhs_pressure(){};

void SOR(){};

void compute_rit(float& rit) {
    rit *= 2;
};

void compute_uv(){};


int main() {
    // Fluid
    float Re = 100;
    float u0 = 0.5f;
    float v0 = 0.0f;
    float p0 = 0.0f;

    // time
    float t = 0;
    float t_max = 100;
    float dt;
    float tau = 0.5;
    int n = 0;

    // grid
    int dims = 2;
    int size_x = 24;
    int size_y = 32;
    float dx = 1.0f;
    float dy = 1.0f;

    int imax = size_x / dx;
    int jmax = size_y / dy;


    // pressure
    int it = 0;
    int it_max = 10;
    float rit = 1.0f;
    float eps = 0.01f;
    float norm_p0 = 0.5f;


    // Initialize the fields
    // TODO: Investigate if order of the dim / grid or if grouping the fields together changes the perf.
    Tensor U({dims, imax+2, jmax+2}, u0);
    Tensor V({dims, imax+2, jmax+2}, v0);
    Tensor P({dims, imax+2, jmax+2}, p0);
    Tensor F({dims, imax+2, jmax+2});
    Tensor G({dims, imax+2, jmax+2});
    Tensor RHS({dims, imax+2, jmax+2});
    Tensor OBS({imax+2, jmax+2});

    // TODO: Remove this hardcoded boundary
    // Set up and down part of box
    enum bd {
        NO,             // Not Boundary
        YES,            // Boundary value
        UP,             // The boundary is above
        DOWN,           // The boundary is below
        LEFT,           // The boundary is to the left
        RIGHT           // The boundary is to the right
    };

    for (int i=0; i!=OBS.shape(0); ++i) {
        OBS({i, 0}) = bd::YES;
        OBS({i, 1}) = bd::DOWN;
        OBS({i, OBS.shape(1)-1}) = bd::YES;
        OBS({i, OBS.shape(1)-2}) = bd::UP;
    }
    for (int j=1; j!=OBS.shape(1)-1; ++j) {
        OBS({0, j}) = bd::YES;
        OBS({1, j}) = bd::LEFT;
        OBS({OBS.shape(0)-1, j}) = bd::YES;
        OBS({OBS.shape(0)-2, j}) = bd::RIGHT;
    }

    std::cout << OBS << std::endl;

    while (t < t_max) {
        dt = adaptive_time_step_size(U, V, dx, dy, Re, tau, dt);

        set_boundary_values(U, V, OBS);
        compute_FG();
        compute_rhs_pressure();
        it = 0;

        while (it<it_max && rit > eps * norm_p0) {
            SOR();
            compute_rit(rit);
            ++it;
        }

        compute_uv();

        t += dt;
        n += 1;
    }

    return 0;
}
