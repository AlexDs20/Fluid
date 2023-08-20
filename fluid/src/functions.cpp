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
    // TODO: Handle gx gy!!
    // F = u + dt * (1./Re * (dudxdx + dudydy) - duudx - duvdy + gx);
    // G = v + dt * (1./Re * (dvdxdx + dvdydy) - duvdx - dvvdy + gy);

    const static int imax = F.imax();
    const static int jmax = F.jmax();

    const float Reinv = 1.0f/Re;
    const static float dxinv = 1.0f/dx;
    const static float dyinv = 1.0f/dy;
    const static float dxinv2 = dxinv*dxinv;
    const static float dyinv2 = dyinv*dyinv;
    const static float dxinv4 = dxinv * 0.25;
    const static float gammadxinv4 = gamma * dxinv4;
    const static float dyinv4 = dyinv * 0.25;
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

    float one, two, three, four;

    for (int i=1; i!=imax+1; ++i) {
        for (int j=1; j!=jmax+1; ++j) {
            uijm = U({i,   j-1});
            uij  = U({i,   j});
            uijp = U({i,   j+1});
            uimj = U({i-1, j});
            uimjp= U({i-1, j+1});
            uipj = U({i+1, j});

            vijm = V({i,   j-1});
            vij  = V({i,   j});
            vijp = V({i,   j+1});
            vimj = V({i-1, j});
            vipjm= V({i+1, j-1});
            vipj = V({i+1, j});

            if (i<imax)
            {
                dudxdx = (uipj - 2*uij + uimj) * dxinv2;
                dudydy = (uijp - 2*uij + uijm) * dyinv2;

                one   = uij + uipj;
                three = uimj + uij;
                duudx =        dxinv4 * (      one  * one -      three  * three )       \
                        + gammadxinv4 * ( fabs(one) * one - fabs(three) * three );

                one   = vij + vipj;
                two   = uij + uijp;
                three = vijm + vipjm;
                four  = uijm + uij;
                duvdy =        dyinv4 * (      one  * two -       three  * four )      \
                        + gammadyinv4 * ( fabs(one) * two -  fabs(three) * four );

                // TODO:
                float gx = 0.0f;
                F({i, j}) = uij + dt * (Reinv * (dudxdx + dudydy) - duudx - duvdy + gx);
            }

            if (j<jmax)
            {
                dvdxdx = (vipj - 2 * vij + vimj) * dxinv2;
                dvdydy = (vijp - 2 * vij + vijm) * dyinv2;

                one   = uij + uijp;
                two   = vij + vipj;
                three = uimj + uimjp;
                four  = vimj + vij;
                duvdx =        dxinv4 * (      one  * two -      three  * four)        \
                        + gammadxinv4 * ( fabs(one) * two - fabs(three) * four );


                one   = vij + vijp;
                three = vijm + vij;
                dvvdy =        dyinv4 * (      one  * one -      three  * three )      \
                        + gammadyinv4 * ( fabs(one) * one - fabs(three) * three ) ;

                // TODO:
                float gy = -10.0f;
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
    static Tensor Res = P;
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

    for (int i=1; i!=imax+1; ++i) {
        // TODO: Handle boundary conditions
        if (i==1)
            eiW = 0;
        else
            eiW = 1;
        if (i==imax)
            eiE = 0;
        else
            eiE = 1;

        for (int j=1; j!=jmax+1; ++j) {
            // TODO: Handle boundary conditions
            if (j==1)
                ejS = 0;
            else
                ejS = 1;
            if (j==jmax)
                ejN = 0;
            else
                ejN = 1;

            pit_imj = pit({i-1, j});
            pit_ipj = pit({i+1, j});
            pit_ijm = pit({i, j-1});
            pit_ij  = pit({i, j});
            pit_ijp = pit({i, j+1});

            P({i, j}) = (1-omega) * pit_ij         \
                + omega / ( dxinv2 * (eiE+eiW) + dyinv2 * (ejN+ejS) )       \
                * ( dxinv2 * ( eiE * pit_ipj + eiW * P({i-1, j}) ) + dyinv2 * (ejN * pit_ijp + ejS * P({i, j-1})) - RHS({i, j}) );

            // Res({i, j})
            rit_tmp = dxinv2 * ( eiE * (pit_ipj - pit_ij) - eiW * (pit_ij - pit_imj) )          \
                + dyinv2 * ( ejN * (pit_ijp - pit_ij) - ejS * (pit_ij - pit_ijm) )                  \
                - RHS({i, j});
            rit = std::max(fabs(rit_tmp), rit);
        }
    }

};

void compute_rit(float& rit) {
    rit *= 2;
};

void compute_uv(){};


