#include <iostream>
#include <cmath>
#include <omp.h>
#include "Physics/calculate.hpp"


void set_constant_flags(Domain& domain, int imax, int jmax) {
    // Set the flags that do not change over time
    // Boundaries and sphere in the middle
    // Sphere info
    const int px = imax / 5;
    const int py = jmax / 2;
    const float radius = std::min(imax, jmax) / 10.0f;

    for (int j=0; j<jmax+2; ++j) {
        for (int i=0; i<imax+2; ++i) {
            // Above and Below
            if (i==0 || i==imax+1) {
                domain(i, j).obstacle = true;
            } else if (j==0 || j==jmax+1) {      // left / right
                domain(i, j).obstacle = true;
            } else if ( (px-i)*(px-i) + (py-j)*(py-j) < radius*radius) {    // Mark the sphere
                domain(i, j).obstacle = true;
            } else {
                domain(i, j).obstacle = false;
            }
        }
    }

    // Set where the water is in the obstacle cells
    for (int j=0; j<jmax+2; ++j) {
        for (int i=0; i<imax+2; ++i) {
            if (domain(i, j).obstacle) {
                if (i<imax+1 & domain(i+1, j).obstacle == false){
                    domain(i, j).E = true;
                }
                if (i>0 & domain(i-1, j).obstacle == false){
                    domain(i, j).W = true;
                }
                if (j<jmax+1 & domain(i, j+1).obstacle == false){
                    domain(i, j).N = true;
                }
                if (j>0 & domain(i, j-1).obstacle == false){
                    domain(i, j).S = true;
                }
            }
        }
    }
};

void set_boundary_values(Tensor& U, Tensor& V, const Domain& domain, int imax, int jmax) {
    for (int j=0; j<jmax+2; ++j) {
        for (int i=0; i<imax+2; ++i) {
            if (domain(i, j).obstacle) {
                // TODO: check boundary type
                //       Start with no slip
                // U
                if (domain(i, j).W) {
                    U(i-1, j) = 0;
                } else {
                    if (domain(i, j).N) {
                        U(i-1, j) = -U(i-1, j+1);
                    }
                    if (domain(i, j).S) {
                        U(i-1, j) = -U(i-1, j-1);
                    }
                }
                if (domain(i, j).E) {
                    U(i, j) = 0;
                } else {
                    if (domain(i, j).N) {
                        U(i, j) = -U(i, j+1);
                    }
                    if (domain(i, j).S) {
                        U(i, j) = -U(i, j-1);
                    }
                }

                // V
                if (domain(i, j).N) {
                    V(i, j) = 0;
                } else {
                    if (domain(i, j).W) {
                        V(i, j) = -V(i-1, j);
                    }
                    if (domain(i, j).E) {
                        V(i, j) = -V(i+1, j);
                    }
                }
                if (domain(i, j).S) {
                    V(i, j-1) = 0;
                } else {
                    if (domain(i, j).W) {
                        V(i, j-1) = -V(i-1, j-1);
                    }
                    if (domain(i, j).E) {
                        V(i, j-1) = -V(i+1, j-1);
                    }
                }
            }
        }
    }
};

void set_specific_boundary_values(Tensor& U, Tensor& V, int imax, int jmax) {
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
    const int width = 4;
    for (int j=(jmax-width)/2; j!=(jmax+width)/2;++j)
    {
        U(0, j) = u;
    }
};

float adaptive_time_step_size( const Tensor& U, const Tensor& V, float dt, const Constants& c) {
    const float dx     = c.dx;
    const float dy     = c.dy;
    const float tau    = c.tau;
    const int imax     = c.imax;
    const int jmax     = c.jmax;
    const float Re_dt  = 0.5f * c.Re / (c.dxdxinv + c.dydyinv);

    if ( (tau <= 0) || (tau > 1) ) {
        return dt;
    }

    float absumax=0, absvmax=0;
    for (int j=0; j<jmax+2; ++j) {
        for (int i=0; i<imax+2; ++i){
            absumax = std::max(absumax, std::fabs(U(i, j)));
            absvmax = std::max(absvmax, std::fabs(V(i, j)));
        }
    }

    const float dt_v_max = std::min(dx / absumax, dy / absvmax);

    float ret = std::min(Re_dt, dt_v_max);
    ret = std::min(dt, tau*ret);

    return ret;
};

void compute_FG(Tensor& F, Tensor& G, const Tensor& U, const Tensor& V, const Domain& domain, float dt, float gx, float gy, const Constants& c) {
    // F = u + dt * (1./Re * (dudxdx + dudydy) - duudx - duvdy + gx);       // i=1..imax-1 j=1..jmax
    // G = v + dt * (1./Re * (dvdxdx + dvdydy) - duvdx - dvvdy + gy);       // i=1..imax   j=1..jmax-1

    // 100 math ops per loop
    // 100 * (iamx+2)*(jamx+2) = 100 * 256 * 128 = 3_276_800 math ops
    // cpu:
    // 2.7Ghz -> 2_700_000_000 clocks / seconds
    //  simd -> * 8
    //  core -> * 4         (*8 if consider logical cores)

    const float Re          = c.Re;
    const float dx          = c.dx;
    const float dy          = c.dy;
    const int imax          = c.imax;
    const int jmax          = c.jmax;
    const float Reinv       = 1.0f/Re;

    const float dxinv       = c.dxinv;
    const float dyinv       = c.dyinv;

    const float dxinv2      = c.dxdxinv;
    const float dyinv2      = c.dydyinv;

    const float dxinv4      = c.dxinv4;
    const float gammadxinv4 = c.gammadxinv4;

    const float dyinv4      = c.dyinv4;
    const float gammadyinv4 = c.gammadyinv4;

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

    // U
    float uijm;
    float uij;
    float uijp;

    float uimj;
    float uimjp;
    float uipj;

    // V
    float vijm;
    float vij;
    float vijp;

    float vimj;
    float vipjm;
    float vipj;

    float one, two, three, four, five, six;

    for (int j=0; j!=jmax+2; ++j) {
        for (int i=0; i!=imax+2; ++i) {
            if (domain(i, j).obstacle == false) {
                uimj = U(i-1, j);
                uimjp= U(i-1, j+1);
                uijm = U(i,   j-1);
                uij  = U(i,   j);
                uijp = U(i,   j+1);
                uipj = U(i+1, j);

                vimj = V(i-1, j);
                vijm = V(i,   j-1);
                vij  = V(i,   j);
                vijp = V(i,   j+1);
                vipjm= V(i+1, j-1);
                vipj = V(i+1, j);

                if (i<imax+1 & domain(i+1, j).obstacle == false)
                {
                    dudxdx = (uipj - 2.0f*uij + uimj) * dxinv2;
                    dudydy = (uijp - 2.0f*uij + uijm) * dyinv2;

                    one   = uij  + uipj;
                    two   = uimj + uij;
                    three = uij  - uipj;
                    four  = uimj - uij;
                    duudx =        dxinv4 * (      one  * one   -      two  * two )       \
                            + gammadxinv4 * ( std::fabs(one) * three - std::fabs(two) * four );
                    // 23

                    one   = vij + vipj;
                    two   = uij + uijp;
                    three = vijm + vipjm;
                    four  = uijm + uij;
                    five = uij - uijp;
                    six =  uijm - uij;
                    duvdy =        dyinv4 * (      one  * two -        three  * four )      \
                            + gammadyinv4 * ( std::fabs(one) * five -  std::fabs(three) * six );
                    // 17

                    F(i, j) = uij + dt * (Reinv * (dudxdx + dudydy) - duudx - duvdy + gx);
                    // 7
                }

                if (j<jmax+1 & domain(i, j+1).obstacle == false)
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
                            + gammadxinv4 * ( std::fabs(one) * five - std::fabs(three) * six );

                    one   = vij  + vijp;
                    three = vijm + vij;
                    five  = vij  - vijp;
                    six   = vijm - vij;
                    dvvdy =        dyinv4 * (      one  * one  -      three  * three )      \
                            + gammadyinv4 * ( std::fabs(one) * five - std::fabs(three) * six ) ;

                    G(i, j) = vij + dt * (Reinv * (dvdxdx + dvdydy) - duvdx - dvvdy + gy );
                    // 48?
                }
            } else {
                if (domain(i, j).N)
                    G(i, j) = V(i, j);
                else if (domain(i, j).S)
                    G(i, j-1) = V(i, j-1);

                if (domain(i, j).W)
                    F(i-1, j) = U(i-1, j);
                else if (domain(i, j).E)
                    F(i, j) = U(i, j);
            }
        }
    }
};

void compute_rhs_pressure(Tensor& RHS, const Tensor& F, const Tensor& G, const Domain& domain, float dt, const Constants& c) {
    const int imax    = c.imax;
    const int jmax    = c.jmax;
    const float dxinv = c.dxinv;
    const float dyinv = c.dyinv;

    const float dtinv = 1.0f / dt;

    for (int j=0; j<jmax+2; ++j) {
        for (int i=0; i<imax+2; ++i) {
            if (domain(i, j).obstacle==false) {
                RHS(i, j) = dtinv * ( dxinv * ( F(i, j) - F(i-1, j) ) + dyinv * ( G(i, j) - G(i, j-1) ) );
            }
        }
    }
};

void SOR(Tensor& P, const Tensor& RHS, const Domain& domain, float& rit, const Constants& c) {
    const int imax     = c.imax;
    const int jmax     = c.jmax;
    const float dxinv2 = c.dxdxinv;
    const float dyinv2 = c.dydyinv;
    const float omega  = c.omega;
    const float coeff  = c.coeff;

    Tensor Pt({jmax+2, imax+2}, 0.0f);


    float rit_tmp;

    rit = 0.0f;

    float tmp_p;
    int edges;

    // Set pressure boundary conditions
    for (int j=0; j!=jmax+2; ++j) {
        for (int i=0; i!=imax+2; ++i) {
            if (domain(i, j).obstacle) {
                tmp_p = 0.0f;
                edges = 0;
                if (domain(i, j).N) {
                    tmp_p += P(i, j+1);
                    ++edges;
                }
                if (domain(i, j).S) {
                    tmp_p += P(i, j-1);
                    ++edges;
                }
                if (domain(i, j).E) {
                    tmp_p += P(i+1, j);
                    ++edges;
                }
                if (domain(i, j).W) {
                    tmp_p += P(i-1, j);
                    ++edges;
                }
                P(i, j) = tmp_p / (float)edges;
            }
        }
    }


    // Update pressure
    for (int j=1; j!=jmax+1; ++j) {
        for (int i=1; i!=imax+1; ++i) {
            if (domain(i, j).obstacle == false) {
                P(i, j) = (1.0f - omega) * P(i, j)                    \
                    + coeff                                           \
                    * ( dxinv2 * ( P(i+1, j) + P(i-1, j) )            \
                      + dyinv2 * ( P(i, j+1) + P(i, j-1) )            \
                      - RHS(i, j)                                     \
                      );
            }
        }
    }

    // Compute residual
    for (int j=1; j!=jmax+1; ++j) {
        for (int i=1; i!=imax+1; ++i) {
            if (domain(i, j).obstacle == false) {
                rit_tmp = dxinv2 * ( P(i-1, j) - 2.0f*P(i, j) + P(i+1, j) ) \
                        + dyinv2 * ( P(i, j-1) - 2.0f*P(i, j) + P(i, j+1) ) \
                        - RHS(i, j);
                rit = std::max(std::fabs(rit_tmp), rit);
            }
        }
    }
};


void compute_uv(Tensor& U, Tensor& V, const Tensor& F, const Tensor& G, const Tensor& P, const Domain& domain, float dt, const Constants& c) {
    const int imax    = c.imax;
    const int jmax    = c.jmax;
    const float dxinv = c.dxinv;
    const float dyinv = c.dyinv;

    const float dtdx = dt * dxinv;
    const float dtdy = dt * dyinv;

    float pipj;
    float pij;
    float pijp;

    for (int j=1; j!=jmax+1; ++j) {
        for (int i=1; i!=imax+1; ++i) {
            pij  = P(i, j);
            pijp = P(i, j+1);
            pipj = P(i+1, j);

            if ( (domain(i, j).obstacle==false) && domain(i+1, j).obstacle==false )
                U(i, j) = F(i, j) - dtdx * (pipj - pij);
            if ( (domain(i, j).obstacle==false) && domain(i, j+1).obstacle==false )
                V(i, j) = G(i, j) - dtdy * (pijp - pij);
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
