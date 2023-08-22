#include <iostream>
#include <math.h>
#include "functions.hpp"

float adaptive_time_step_size( const Tensor& U, const Tensor& V, float dx, float dy, float Re, float tau, float dt) {
    const static float dxinv2 = 1.0f / (dx*dx);
    const static float dyinv2 = 1.0f / (dy*dy);
    const static float Re_dt = ( Re / (2.0f * (dxinv2 + dyinv2)) );

    float ret = std::min(
            dt,
            Re_dt
            );

    float dt_v_max = std::min(dx / U.amax(), dy / V.amax());

    return tau * std::min(ret, dt_v_max);
};

void set_boundary_values(Tensor& U, Tensor& V){
    // Currently hard-code up/down left/right
    int imax = U.imax();
    int jmax = U.jmax();

    // Top/Bottom: NO SLIP
    for (int i=1; i!=imax+1; ++i) {
        V({i, 0}) = 0.0f;
        V({i, jmax}) = 0.0f;
        U({i, 0}) = -U({i, 1});
        U({i, jmax+1}) = -U({i, jmax});
    }

    // Left: input flow, Right: output
    float u = 0.7;
    float v = 0.0;
    for (int j=1; j!=jmax+1; ++j) {
        // U({0, j}) = u;
        // V({0, j}) = v;
        U({0, j}) = 0;
        V({0, j}) = -V({1, j});

        U({imax, j}) = U({imax-1, j});
        V({imax+1, j}) = V({imax, j});
    }

    for (int j=(jmax+1)/2-2; j!=(jmax+1)/2+3;++j)
        U({0, j}) = u;
};

void set_object_boundary_values(Tensor& U, Tensor& V, Tensor& OBS) {
    // TODO
};

void compute_FG(Tensor& F, Tensor& G, const Tensor& U, const Tensor& V, float dt, float Re, float dx, float dy, float gamma) {
    // TODO: Handle gx gy!!
    // F = u + dt * (1./Re * (dudxdx + dudydy) - duudx - duvdy + gx);       // i=1..imax-1 j=1..jmax
    // G = v + dt * (1./Re * (dvdxdx + dvdydy) - duvdx - dvvdy + gy);       // i=1..imax   j=1..jmax-1

    const static int imax = F.imax();
    const static int jmax = F.jmax();

    const static float Reinv = 1.0f/Re;

    const static float dxinv = 1.0f/dx;
    const static float dyinv = 1.0f/dy;

    const static float dxinv2 = dxinv*dxinv;
    const static float dyinv2 = dyinv*dyinv;

    const static float dxinv4 = dxinv * 0.25f;
    const static float gammadxinv4 = gamma * dxinv4;

    const static float dyinv4 = dyinv * 0.25f;
    const static float gammadyinv4 = gamma * dyinv4;

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
    float uimjp;
    float uipj;

    float vijm;
    float vij;
    float vijp;

    float vimj;
    float vipjm;
    float vipj;

    float one, two, three, four, five, six;

    for (int i=1; i!=imax+1; ++i) {
        for (int j=1; j!=jmax+1; ++j) {
            uimj = U({i-1, j});
            uimjp= U({i-1, j+1});
            uijm = U({i,   j-1});
            uij  = U({i,   j});
            uijp = U({i,   j+1});
            uipj = U({i+1, j});

            vimj = V({i-1, j});
            vijm = V({i,   j-1});
            vij  = V({i,   j});
            vijp = V({i,   j+1});
            vipjm= V({i+1, j-1});
            vipj = V({i+1, j});

            //if (i<imax)
            {
                dudxdx = (uipj - 2*uij + uimj) * dxinv2;
                dudydy = (uijp - 2*uij + uijm) * dyinv2;

                one   = uij  + uipj;
                two   = uimj + uij;
                three = uij  - uipj;
                four  = uimj - uij;
                duudx =        dxinv4 * (      one  * one   -      two  * two )       \
                        + gammadxinv4 * ( fabs(one) * three - fabs(two) * four );

                one   = vij + vipj;
                two   = uij + uijp;
                three = vijm + vipjm;
                four  = uijm + uij;
                five = uij - uijp;
                six =  uijm - uij;
                duvdy =        dyinv4 * (      one  * two -        three  * four )      \
                        + gammadyinv4 * ( fabs(one) * five -  fabs(three) * six );

                // TODO:
                float gx = 0.0f;
                F({i, j}) = uij + dt * (Reinv * (dudxdx + dudydy) - duudx - duvdy + gx);
            }

            //if (j<jmax)
            {
                dvdxdx = (vipj - 2 * vij + vimj) * dxinv2;
                dvdydy = (vijp - 2 * vij + vijm) * dyinv2;

                one   = uij  + uijp;
                two   = vij  + vipj;
                three = uimj + uimjp;
                four  = vimj + vij;
                five  = vij  - vipj;
                six   = vimj - vij;
                duvdx =        dxinv4 * (      one  * two  -      three  * four)        \
                        + gammadxinv4 * ( fabs(one) * five - fabs(three) * six );

                one   = vij  + vijp;
                three = vijm + vij;
                five  = vij  - vijp;
                six   = vijm - vij;
                dvvdy =        dyinv4 * (      one  * one  -      three  * three )      \
                        + gammadyinv4 * ( fabs(one) * five - fabs(three) * six ) ;

                // TODO:
                float gy = 0.0f;
                G({i, j}) = vij + dt * (Reinv * (dvdxdx + dvdydy) - duvdx - dvvdy + gy );
            }
        }
    }
};

void compute_rhs_pressure(Tensor& RHS, const Tensor& F, const Tensor& G, float dx, float dy, float dt) {
    const static float dxinv = 1.0f/dx;
    const static float dyinv = 1.0f/dy;
    const static float imax = F.imax();
    const static float jmax = F.jmax();

    float dtinv = 1.0f / dt;

    for (int i=1; i!=imax+1; ++i)
        for (int j=1; j!=jmax+1; ++j)
            RHS({i, j}) = dtinv * ( dxinv * (F({i, j}) - F({i-1, j})) + dyinv * (G({i, j}) - G({i, j-1})) );
};

void SOR(Tensor& P, float& rit, const Tensor& RHS, float omega, float dx, float dy) {
    // static Tensor Res = P;
    static Tensor pit = P;
    pit = P;

    const static int imax = P.imax();
    const static int jmax = P.jmax();
    const static float dxinv2 = (1.0f / dx) * (1.0f / dx);
    const static float dyinv2 = (1.0f / dy) * (1.0f / dy);

    float rit_tmp;

    float pit_imj;
    float pit_ij;
    float pit_ijm;
    float pit_ipj;
    float pit_ijp;

    float eiE;
    float eiW;
    float ejN;
    float ejS;

    rit = 0.0f;

    for (int i=1; i!=imax+1; ++i) {
        // TODO: Handle boundary conditions
        eiW = 1;
        if (i==1)
            eiW = 0;

        eiE = 1;
        if (i==imax)
            eiE = 0;

        for (int j=1; j!=jmax+1; ++j) {
            // TODO: Handle boundary conditions
            ejS = 1;
            if (j==1)
                ejS = 0;

            ejN = 1;
            if (j==jmax)
                ejN = 0;

            pit_imj = pit({i-1, j});
            pit_ijm = pit({i, j-1});
            pit_ij  = pit({i, j});
            pit_ijp = pit({i, j+1});
            pit_ipj = pit({i+1, j});

            P({i, j}) = (1-omega) * pit_ij         \
                + omega / ( dxinv2 * (eiE+eiW) + dyinv2 * (ejN+ejS) )       \
                * ( dxinv2 * ( eiE * pit_ipj + eiW * P({i-1, j}) ) + dyinv2 * (ejN * pit_ijp + ejS * P({i, j-1})) - RHS({i, j}) );

            // Res({i, j})
            rit_tmp = dxinv2 * ( eiE * (pit_ipj - pit_ij) - eiW * (pit_ij - pit_imj) )              \
                    + dyinv2 * ( ejN * (pit_ijp - pit_ij) - ejS * (pit_ij - pit_ijm) )              \
                    - RHS({i, j});
            rit = std::max(fabs(rit_tmp), rit);
        }
    }

};

void compute_uv(Tensor& U, Tensor& V, const Tensor& F, const Tensor& G, const Tensor& P, float dx, float dy, float dt) {
    const static float dxinv = 1.0f / dx;
    const static float dyinv = 1.0f / dy;

    const static float imax = U.imax();
    const static float jmax = U.jmax();

    float dtdx = dt * dxinv;
    float dtdy = dt * dyinv;

    float pipj;
    float pij;
    float pijp;

    for (int i=1; i!=imax+1; ++i) {
        for (int j=1; j!=jmax+1; ++j) {
            pij  = P({i, j});
            pijp = P({i, j+1});
            pipj = P({i+1, j});

            if (i!=imax)
                U({i, j}) = F({i, j}) - dtdx * (pipj - pij);
            if (j!=jmax)
                V({i, j}) = G({i, j}) - dtdy * (pijp - pij);
        }
    }
};


