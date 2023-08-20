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
        V({i, jmax}) = 0.0f;
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

void compute_FG(Tensor& F, Tensor& G, const Tensor& U, const Tensor& V, float dt, float Re, float dx, float dy, float gamma) {
    int imax = F.imax();
    int jmax = F.jmax();

    float Reinv = 1.0f/Re;
    float dxinv = 1.0f/dx;
    float dyinv = 1.0f/dy;

    // F
    float dudxdx;
    float dudydy;
    float duudx;
    float duvdy;

    // G
    float dvdxdx;
    float dvdydy;
    float duvdx;
    float dvvdy;

    // Good for cache?
    float uijm;
    float uij;
    float uijp;

    float uimj;
    float uipj;

    float vijm;
    float vij;
    float vijp;

    float vimj;
    float vipjm;
    float vipj;

    float one, two, three, four;

    for (int i=1; i!=imax; ++i) {
        for (int j=1; j!=jmax; ++j) {
            uijm = U({i,   j-1});
            uij  = U({i,   j});
            uijp = U({i,   j+1});
            uimj = U({i-1, j});
            uipj = U({i+1, j});

            vijm = V({i,   j-1});
            vij  = V({i,   j});
            vijp = V({i,   j+1});
            vimj = V({i-1, j});
            vipjm= V({i+1, j-1});
            vipj = V({i+1, j});

            dudxdx = (uipj - 2*uij + uimj) * dxinv;
            dudydy = (uijp - 2*uij + uijm) * dyinv;

            one   = uij + uipj;
            three = uimj + uij;
            duudx = dxinv * 0.25           * (      one  * one -      three  * three )       \
                    + gamma * dxinv * 0.25 * ( fabs(one) * one - fabs(three) * three );

            one   = vij + vipj;
            two   = uij + uijp;
            three = vijm + vipjm;
            four  = uijm + uij;
            duvdy = dyinv * 0.25           * (      one  * two -       three  * four )      \
                    + gamma * dyinv * 0.25 * ( fabs(one) * two -  fabs(three) * four );

            F({i, j}) = uij + dt * (Reinv * (dudxdx + dudydy) - duudx - duvdy );
        }
    }

    // F = u + dt * (1./Re * (dudxdx + dudydy) - duudx - duvdy + gx);
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


