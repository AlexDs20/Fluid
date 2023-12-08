#include <iostream>
#include <cmath>
#include <omp.h>
#include "Math/Matrix.hpp"
#include "Physics/calculate.hpp"
#include "Physics/define.hpp"
#include "Math/simd.h"


void set_obstacle_flags(Matrixi& domain, int imax, int jmax) {
    const int k = imax+2;
    const int l = jmax+2;
    const wide_float px((float)k/5);
    const wide_float py((float)l/2);
    const wide_float r((float) l*l/100);

    const wide_int zero(0);
    const wide_int obs(OBSTACLE);
    const wide_int lm1(l-1);
    const wide_int km1(k-1);

    for (int j=0; j<l; j++) {
        // for j
        wide_int wide_j(j);
        wide_int wide_mask_j = (wide_j == zero) | (wide_j == lm1);

        // Sphere
        wide_float wide_sphere_j = Square(py-wide_j);

        for (int i=0; i+7<k; i+=8) {
            // for i
            wide_int wide_i = WideIndex(i);
            wide_int wide_mask_i = (wide_i == 0) | (wide_i == km1);

            // Combine i and j masks
            wide_int wide_mask = wide_mask_j | wide_mask_i;

            // Sphere
            wide_int wide_sphere = (wide_sphere_j + Square(px - wide_i)) < r;

            // Combine sphere mask and ij mask
            wide_mask |= wide_sphere;
            wide_int assign_value(&domain(i, j));
            ConditionalAssign(&assign_value, wide_mask, obs);

            StoreWideInt(&domain(i, j), assign_value);
        }
    }
};

void set_fluid_flags(Matrixi& domain, int imax, int jmax) {
    const wide_int obs(OBSTACLE);
    const wide_int north(N);
    const wide_int south(S);
    const wide_int west(W);
    const wide_int east(E);

    const wide_int ip1(imax+1);
    const wide_int jp1(jmax+1);
    const wide_int zero(0);

    // Set where the water is in the obstacle cells
    for (int j=0; j<jmax+2; ++j) {
        wide_int wide_j(j);

        for (int i=0; i+7<imax+2; i+=8) {
            wide_int wide_i = WideIndex(i);

            wide_int domij(&domain(i, j));

            // TODO(alex): HOW TO HANDLE THIS?! Go Over the domain size!?
            if (endi) {
                wide_int ip1mask = LoadInts(-1,-1,-1,-1,-1,-1,-1,0);
                wide_int domipj = LoadMaskedPackedWideInt(&domain(i+1, j), ip1mask);
            } else {
                wide_int domipj(&domain(i+1, j));
            }
            wide_int domimj(&domain(i-1, j));

            wide_int domijp(&domain(i, j+1));
            wide_int domijm(&domain(i, j-1));

            wide_int mask = (domij == obs);
            // North
            wide_int mask_neighbour = ((wide_j<jp1) & ((domijp & obs) == zero));
            domij |= (north & (mask & mask_neighbour));
            // South
            mask_neighbour = ((wide_j>0) & ((domijm & obs) == zero));
            domij |= (south & (mask & mask_neighbour));
            // West
            mask_neighbour = ((wide_i > zero) & ((domimj & obs) == zero));
            domij |= (west & (mask & mask_neighbour));
            // East
            mask_neighbour = ((wide_i < ip1) & ((domipj & obs) == zero));
            domij |= (east & (mask & mask_neighbour));

            StoreWideInt(&domain(i, j), domij);
        }
    }
};

void set_constant_flags(Matrixi& domain, int imax, int jmax) {
    set_obstacle_flags(domain, imax, jmax);
    set_fluid_flags(domain, imax, jmax);
};


void set_boundary_values(Matrix& U, Matrix& V, const Matrixi& domain, int imax, int jmax) {
    for (int j=0; j<jmax+2; ++j) {
        for (int i=0; i<imax+2; ++i) {
            if ((domain(i, j) & OBSTACLE) == 1) {
                // TODO: check boundary type
                //       Start with no slip
                // U
                if ((domain(i, j) & W)) {
                    U(i-1, j) = 0;
                } else {
                    if ((domain(i, j) & N)) {
                        U(i-1, j) = -U(i-1, j+1);
                    }
                    if ((domain(i, j) & S)) {
                        U(i-1, j) = -U(i-1, j-1);
                    }
                }
                if ((domain(i, j) & E)) {
                    U(i, j) = 0;
                } else {
                    if ((domain(i, j) & N)) {
                        U(i, j) = -U(i, j+1);
                    }
                    if ((domain(i, j) & S)) {
                        U(i, j) = -U(i, j-1);
                    }
                }

                // V
                if ((domain(i, j) & N)) {
                    V(i, j) = 0;
                } else {
                    if ((domain(i, j) & W)) {
                        V(i, j) = -V(i-1, j);
                    }
                    if ((domain(i, j) & E)) {
                        V(i, j) = -V(i+1, j);
                    }
                }
                if ((domain(i, j) & S)) {
                    V(i, j-1) = 0;
                } else {
                    if ((domain(i, j) & W)) {
                        V(i, j-1) = -V(i-1, j-1);
                    }
                    if ((domain(i, j) & E)) {
                        V(i, j-1) = -V(i+1, j-1);
                    }
                }
            }
        }
    }
};

void set_specific_boundary_values(Matrix& U, Matrix& V, int imax, int jmax) {
    // Outflow
    for (int j=1; j!=jmax+1; ++j) {
        U(imax, j) = U(imax-1, j);
        V(imax+1, j) = V(imax, j);
    }
    for (int i=1; i!=imax+1; ++i) {
        U(i, jmax+1) = U(i, jmax);
        V(i, jmax) = V(i, jmax-1);

        V(i, 0) = V(i, 1);
        U(i, 0) = U(i, 1);
    }

    // Left: input flow
    const float u = 0.8;
    const int width = 32;
    for (int j=(jmax-width)/2; j!=(jmax+width)/2;++j)
    {
        U(0, j) = u;
    }
};


float adaptive_time_step_size(Matrix& U, Matrix& V, float dt, const Constants& c) {
    const float dx     = c.dx;
    const float dy     = c.dy;
    const float tau    = c.tau;
    const int imax     = c.imax;
    const int jmax     = c.jmax;
    const float Re_dt  = 0.5f * c.Re / (c.dxdxinv + c.dydyinv);

    if ( (tau <= 0) || (tau > 1) ) {
        return dt;
    }

    wide_float absu(0.0f);
    wide_float absv(0.0f);

    for (int j=0; j<jmax+2; ++j) {
        for (int i=0; i<imax+2; i+=8){
            wide_float uij(&U(i, j));
            wide_float vij(&V(i, j));

            wide_float local_absu = Abs(uij);
            wide_float local_absv = Abs(vij);

            absu = Max(local_absu, absu);
            absv = Max(local_absv, absv);
        }
    }

    // combine U and V maxes
    float absumax = HorizontalMax(absu);
    float absvmax = HorizontalMax(absv);

    const float dt_v_max = std::min(dx / absumax, dy / absvmax);

    float ret = std::min(Re_dt, dt_v_max);
    ret = std::min(dt, tau*ret);

    return ret;
};


void compute_FG(Matrix& F, Matrix& G, Matrix& U, Matrix& V, Matrixi& domain, float dt, float gx, float gy, const Constants& c) {
    // F = u + dt * (1./Re * (dudxdx + dudydy) - duudx - duvdy + gx);       // i=1..imax-1 j=1..jmax
    // G = v + dt * (1./Re * (dvdxdx + dvdydy) - duvdx - dvvdy + gy);       // i=1..imax   j=1..jmax-1

    // 100 math ops per loop
    // 100 * (iamx+2)*(jamx+2) = 100 * 256 * 128 = 3_276_800 math ops
    // cpu:
    // 2.7Ghz -> 2_700_000_000 clocks / seconds
    //  simd -> * 8
    //  core -> * 4         (*8 if consider logical cores)

    const float Re          = c.Re;
    const int imax          = c.imax;
    const int jmax          = c.jmax;
    const float Reinv       = 1.0f/Re;

    const float dxinv2      = c.dxdxinv;
    const float dyinv2      = c.dydyinv;

    const float dxinv4      = c.dxinv4;
    const float gammadxinv4 = c.gammadxinv4;

    const float dyinv4      = c.dyinv4;
    const float gammadyinv4 = c.gammadyinv4;

    // F
    wide_float dudxdx;
    wide_float dudydy;
    wide_float duudx;
    wide_float duvdy;

    // G
    wide_float dvdxdx;
    wide_float dvdydy;
    wide_float duvdx;
    wide_float dvvdy;

    wide_float one, two, three, four, five, six;

    wide_int imp1(imax+1);
    wide_int jmp(jmax+1);

    for (int j=0; j!=jmax+2; ++j) {
        wide_int wide_j(j);

        for (int i=0; i+7<imax+2; i+=8) {
            wide_float uimj(&U(i-1, j));
            wide_float uimjp(&U(i-1, j+1));
            wide_float uijm(&U(i,   j-1));
            wide_float uij (&U(i,   j));
            wide_float uijp(&U(i,   j+1));
            wide_float uipj(&U(i+1, j));

            wide_float vimj(&V(i-1, j));
            wide_float vijm(&V(i,   j-1));
            wide_float vij (&V(i,   j));
            wide_float vijp(&V(i,   j+1));
            wide_float vipjm(&V(i+1, j-1));
            wide_float vipj(&V(i+1, j));

            wide_int domij(&domain(i, j));
            wide_int domipj(&domain(i+1, j));
            wide_int domijp(&domain(i, j+1));
            wide_int obs(OBSTACLE);

            // if ((domain(i, j) & OBSTACLE) == 0) {
            //      if ((i<imax+1) & ((domain(i+1, j) & OBSTACLE) == 0))    F(i, j) = ;
            //      if ((j<jmax+1) & ((domain(i, j+1) & OBSTACLE) == 0))    G(i, j) = ;
            // } else {
            //     if ((domain(i, j) & E))          F(i, j) = U(i, j);
            //     if ((domain(i, j) & W))     F(i-1, j) = U(i-1, j);
            //
            //     if ((domain(i, j) & N))          G(i, j) = V(i, j);
            //     if ((domain(i, j) & S))     G(i, j-1) = V(i, j-1);
            // }

            wide_float fij(&F(i, j));
            wide_int right_is_fluid_mask = ( (domij & obs) != 0) & ((domij & E)!=0);
            ConditionalAssign(&fij, right_is_fluid_mask, uij);

            wide_int wide_i = WideIndex(i);
            wide_int imp(imax+1);
            wide_int fluid_and_right_is_fluid = ( (domij & obs) == 0 ) & ( (wide_i < imp) & ((domipj & obs) == 0) );

            {
                dudxdx = (uipj - 2.0f*uij + uimj) * dxinv2;
                dudydy = (uijp - 2.0f*uij + uijm) * dyinv2;

                one   = uij  + uipj;
                two   = uimj + uij;
                three = uij  - uipj;
                four  = uimj - uij;
                duudx =        dxinv4 * (      one  * one   -      two  * two )       \
                        + gammadxinv4 * ( Abs(one) * three - Abs(two) * four );
                // 23

                one   = vij + vipj;
                two   = uij + uijp;
                three = vijm + vipjm;
                four  = uijm + uij;
                five  = uij - uijp;
                six   = uijm - uij;
                duvdy = dyinv4 * (      one  * two -        three  * four )      \
                        + gammadyinv4 * ( Abs(one) * five -  Abs(three) * six );
                // 17
                wide_float result = uij + dt * (Reinv * (dudxdx + dudydy) - duudx - duvdy + gx);
                ConditionalAssign(&fij, fluid_and_right_is_fluid, result);
                // 7
            }

            StoreWideFloat(&F(i, j), fij);

            wide_int left_is_fluid_mask = ( (domij & obs) != 0 ) & (domij & W);
            wide_float fimj(&F(i-1, j));
            ConditionalAssign(&fimj, left_is_fluid_mask, uimj);
            StoreWideFloat(&F(i-1, j), fimj);


            wide_float gij(&G(i, j));
            wide_int north_is_fluid_mask = ( (domij & obs) != 0 ) & ((domij & N) != 0);
            ConditionalAssign(&gij, north_is_fluid_mask, vij);

            wide_int current_and_north_is_fluid = ( (domij & obs) == 0 ) & ( (wide_j < jmp) & ((domijp & obs)==0) );

            {
                dvdxdx = (vipj - 2.0f*vij + vimj) * dxinv2;
                dvdydy = (vijp - 2.0f*vij + vijm) * dyinv2;

                one   = uij  + uijp;
                two   = vij  + vipj;
                three = uimj + uimjp;
                four  = vimj + vij;
                five  = vij  - vipj;
                six   = vimj - vij;
                duvdx =        dxinv4 * (      one  * two  -      three  * four)        \
                        + gammadxinv4 * ( Abs(one) * five - Abs(three) * six );

                one   = vij  + vijp;
                three = vijm + vij;
                five  = vij  - vijp;
                six   = vijm - vij;
                dvvdy =        dyinv4 * (      one  * one  -      three  * three )      \
                        + gammadyinv4 * ( Abs(one) * five - Abs(three) * six ) ;

                wide_float result;
                result = vij + dt * (Reinv * (dvdxdx + dvdydy) - duvdx - dvvdy + gy );
                ConditionalAssign(&gij, current_and_north_is_fluid, result);
            }

            StoreWideFloat(&G(i, j), gij);

            wide_int south_is_fluid = ((domij & obs) != 0) & ((domij & S) != 0);
            wide_float gijm(&G(i, j-1));
            ConditionalAssign(&gijm, south_is_fluid, vijm);
            StoreWideFloat(&G(i, j-1), gijm);

        }
    }
};

void compute_rhs_pressure(Matrix& RHS, Matrix& F, Matrix& G, Matrixi& domain, float dt, const Constants& c) {
    const int imax    = c.imax;
    const int jmax    = c.jmax;
    // const float dxinv = c.dxinv;
    // const float dyinv = c.dyinv;
    // const float dtinv = 1.0f / dt;

    wide_float dxinv(c.dxinv);
    wide_float dyinv(c.dyinv);
    wide_float dtinv(1.0f/dt);

    for (int j=0; j<jmax+2; ++j) {
        for (int i=0; i+7<imax+2; i+=8) {
            wide_int domij(&domain(i, j));
            wide_float fij(&F(i, j));
            wide_float fimj(&F(i-1, j));
            wide_float gij(&G(i, j));
            wide_float gijm(&G(i, j-1));

            wide_float result;
            result = dtinv * (dxinv * (fij - fimj) + dyinv * (gij - gijm));
            StoreWideFloat(&RHS(i, j), result);

            // if ((domain(i, j) & OBSTACLE) == 0) {
            //     RHS(i, j) = dtinv * ( dxinv * ( F(i, j) - F(i-1, j) ) + dyinv * ( G(i, j) - G(i, j-1) ) );
            // }
        }
    }
};

void SOR(Matrix& P, const Matrix& RHS, const Matrixi& domain, float& rit, const Constants& c) {
    const int imax     = c.imax;
    const int jmax     = c.jmax;
    const float dxinv2 = c.dxdxinv;
    const float dyinv2 = c.dydyinv;
    const float omega  = c.omega;
    const float coeff  = c.coeff;

    // static Matrix Pt({jmax+2, imax+2}, 0.0f);

    // for (int j=0; j!=jmax+2; ++j)
    //     for (int i=0; i!=imax+2; ++i)
    //         Pt(j, i) = P(i, j);

    float rit_tmp;

    rit = 0.0f;

    float tmp_p;
    int edges;

    // Set pressure boundary conditions
    for (int j=0; j!=jmax+2; ++j) {
        for (int i=0; i!=imax+2; ++i) {
            if (domain(i, j) & OBSTACLE) {
                tmp_p = 0.0f;
                edges = 0;
                if ((domain(i, j) & N)) {
                    // tmp_p += Pt(j+1, i);
                    tmp_p += P(i, j+1);
                    ++edges;
                }
                if ((domain(i, j) & S)) {
                    // tmp_p += Pt(j-1, i);
                    tmp_p += P(i, j-1);
                    ++edges;
                }
                if ((domain(i, j) & E)) {
                    tmp_p += P(i+1, j);
                    ++edges;
                }
                if ((domain(i, j) & W)) {
                    tmp_p += P(i-1, j);
                    ++edges;
                }
                P(i, j) = tmp_p / (float)edges;
                // Pt(j, i) = P(i, j);
            }
        }
    }


    // Update pressure
    for (int j=1; j!=jmax+1; ++j) {
        for (int i=1; i!=imax+1; ++i) {
            if ((domain(i, j) & OBSTACLE) == 0) {
                P(i, j) = (1.0f - omega) * P(i, j)                    \
                    + coeff                                           \
                    * ( dxinv2 * ( P(i-1, j) + P(i+1, j) )            \
                      + dyinv2 * ( P(i, j-1) + P(i, j+1) )            \
                      - RHS(i, j)                                     \
                      );
                      // + dyinv2 * ( Pt(j-1, i) + Pt(j+1, i) )
                // Pt(j, i) = P(i, j);
            }
        }
    }

    // Compute residual
    for (int j=1; j!=jmax+1; ++j) {
        for (int i=1; i!=imax+1; ++i) {
            if ((domain(i, j) & OBSTACLE) == 0) {
                rit_tmp = dxinv2 * ( P(i-1, j) - 2.0f*P(i, j) + P(i+1, j) ) \
                        + dyinv2 * ( P(i, j-1) - 2.0f*P(i, j) + P(i, j+1) ) \
                        - RHS(i, j);
                        // + dyinv2 * ( Pt(j-1, i) - 2.0f*P(i, j) + Pt(j+1, i) )
                rit = std::max(std::fabs(rit_tmp), rit);
            }
        }
    }
};


void compute_uv(Matrix& U, Matrix& V, Matrix& F, Matrix& G, Matrix& P, Matrixi& domain, float dt, const Constants& c) {
    const int imax    = c.imax;
    const int jmax    = c.jmax;
    const float dxinv = c.dxinv;
    const float dyinv = c.dyinv;

    const float dtdx = dt * dxinv;
    const float dtdy = dt * dyinv;

    for (int j=1; j!=jmax+1; ++j) {
        for (int i=1; i+7<=imax+1; i+=8) {
            wide_float pij(&P(i, j));
            wide_float pijp(&P(i, j+1));
            wide_float pipj(&P(i+1, j));

            wide_int dij(&domain(i, j));
            wide_int dipj(&domain(i+1, j));
            wide_int dijp(&domain(i, j+1));
            wide_float fij(&F(i, j));
            wide_float gij(&G(i, j));

            // U
            wide_float result = fij - dtdx * (pipj - pij);
            StoreWideFloat(&U(i, j), result);

            // V
            result = gij - dtdy * (pijp - pij);
            StoreWideFloat(&V(i, j), result);
        }
    }
};

void get_parameters(const std::string& problem, Parameters& params, Constants& constants){
    if (problem == "inflow"){
        // Fluid
        params.Re = 10000;
        params.u0 = 0.0f;
        params.v0 = 0.0f;
        params.p0 = 0.0f;
        params.gx = 0.0f;
        params.gy = 0.0f;

        // time
        params.t_max = 100;
        params.dt_max = 0.04;
        params.tau = 0.3;

        // grid
        params.imax = 256-2;
        params.jmax = 128-2;
        params.xlength = 2.0f;
        params.ylength = 1.0f;

        // pressure
        params.it_max = 20;
        params.eps = 0.001f;
        params.omega = 1.5;

        // Derivative scheme
        params.gamma = 0.9f;
    }

    constants.imax = params.imax;
    constants.jmax = params.jmax;
    constants.dx = params.xlength / params.imax;
    constants.dy = params.ylength / params.jmax;
    constants.dxinv = 1.0f / constants.dx;
    constants.dyinv = 1.0f / constants.dy;
    constants.dxdxinv = constants.dxinv * constants.dxinv;
    constants.dydyinv = constants.dyinv * constants.dyinv;
    constants.dxinv4 = constants.dxinv * 0.25f;
    constants.dyinv4 = constants.dyinv * 0.25f;
    constants.gammadxinv4 = params.gamma * constants.dxinv4;
    constants.gammadyinv4 = params.gamma * constants.dyinv4;
    constants.Re = params.Re;
    constants.tau = params.tau;
    constants.omega = params.omega;
    constants.coeff = params.omega / ( 2.0f * (constants.dxdxinv + constants.dydyinv) );
}
