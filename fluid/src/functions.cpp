#include <iostream>
#include <math.h>
#include "functions.hpp"

float adaptive_time_step_size( const Tensor& U, const Tensor& V, float dx, float dy, float Re, float tau, float dt) {
    float ret = std::min(
            dt,
            (Re / 2.0f) * 1.0f / (1.0f/(dx*dx) + 1.0f/(dy*dy))
            );

    float dt_v_max = std::min(dx / U.amax(), dy / V.amax());

    return tau * std::min(ret, dt_v_max);
};

void set_boundary_values(Tensor& U, Tensor& V){
    // Currently hard-code up/down left/right
    int imax = U.imax();
    int jmax = U.jmax();

    // Top/Bottom: NO SLIP
    for (int i=0; i!=U.shape(0); ++i) {
        V({i, 0}) = 0.0f;
        V({i, jmax-1}) = 0.0f;
        U({i, 0}) = -U({i, 1});
        U({i, jmax}) = -U({i, jmax-1});
    }

    // Left: input flow, Right: output
    float u = 0.2;
    float v = 0.0;
    for (int j=1; j!=U.shape(1)-1; ++j) {
        U({0, j}) = u;
        V({0, j}) = v;

        U({imax, j}) = U({imax-1, j});
        V({imax+1, j}) = V({imax, j});
    }
};

void set_object_boundary_values(Tensor& U, Tensor& V, Tensor& OBS) {
    // TODO
};

void compute_FG(Tensor& F, Tensor& G, const Tensor& U, const Tensor& V, float dt, float Re) {
    float dudxdx = 0.0f;
    float dudydy = 0.0f;
    float duudx = 0.0f;
    float duvdy = 0.0f;

    // F = u + dt * (1./Re * (dudxdx + dudydy) - duudx - duvdy + gx);


    float dvdxdx = 0.0f;
    float dvdydy = 0.0f;
    float duvdx = 0.0f;
    float dvvdy = 0.0f;
    // G = v + dt * (1./Re * (dvdxdx + dvdydy) - duvdx - dvvdy + gy);
};

void compute_rhs_pressure(Tensor& RHS, const Tensor& F, const Tensor& G, float dx, float dy, float dt) {
    int imax = RHS.shape(0)-2;
    int jmax = RHS.shape(1)-2;

    for (int i=1; i!=imax+1; ++i)
        for (int j=1; j!=jmax+1; ++j)
            RHS({i, j}) = (1/dt) * ( ( F({i,j}) - F({i-1,j}) )/dx + ( G({i, j}) - G({i, j-1}) )/dy );
};

void SOR(){
};

void compute_rit(float& rit) {
    rit *= 2;
};

void compute_uv(){};


