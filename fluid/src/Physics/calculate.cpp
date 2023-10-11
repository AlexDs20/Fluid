#include <iostream>
#include <cmath>
#include <omp.h>
#include "Physics/calculate.hpp"


float adaptive_time_step_size( const Tensor& U, const Tensor& V, float dx, float dy, float Re, float tau, float dt, int imax, int jmax) {
    const static float dxinv2 = 1.0f / (dx*dx);
    const static float dyinv2 = 1.0f / (dy*dy);
    const static float Re_dt = 0.5f * Re / (dxinv2 + dyinv2);

    if ( (tau <= 0) || (tau > 1) ) {
        return dt;
    }

    float absumax=0, absvmax=0;
#ifdef _OMP
#pragma omp parallel for reduction(max: absumax, absvmax) default(none) shared(imax, jmax, U, V, std::cout)
#endif
    for (int i=0; i<imax+2; ++i){
        for (int j=0; j<jmax+2; ++j) {
            absumax = std::max(absumax, std::fabs(U(i, j)));
            absvmax = std::max(absvmax, std::fabs(V(i, j)));
        }
    }
#ifdef _OMP
#pragma omp barrier
#endif

    const float dt_v_max = std::min(dx / absumax, dy / absvmax);

    float ret = std::min(Re_dt, dt_v_max);
    ret = std::min(dt, tau*ret);

    return ret;
};

void set_constant_flags(Domain& domain, int imax, int jmax) {
    // Set the flags that do not change over time
    // Boundaries and sphere in the middle
    // Sphere info
    const int px = (imax+2) / 5;
    const int py = (jmax+2) / 2;
    const float radius = std::min(imax, jmax) / 10.0;

#ifdef _OMP
#pragma omp parallel for default(none) shared(domain, imax, jmax, px, py, radius)
#endif
    for (int i=0; i<imax+2; ++i) {
        for (int j=0; j<jmax+2; ++j) {
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
#ifdef _OMP
#pragma omp barrier
#endif

    // Set where the water is in the obstacle cells
#ifdef _OMP
#pragma omp parallel for default(none) shared(domain, imax, jmax)
#endif
    for (int i=0; i<imax+2; ++i) {
        for (int j=0; j<jmax+2; ++j) {
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
#ifdef _OMP
#pragma omp barrier
#endif
};

void set_boundary_values(Tensor& U, Tensor& V, const Domain& domain, int imax, int jmax) {
#ifdef _OMP
#pragma omp parallel for default(none) shared(domain, imax, jmax, U, V)
#endif
    for (int i=0; i<imax+2; ++i) {
        for (int j=0; j<jmax+2; ++j) {
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
#ifdef _OMP
#pragma omp barrier
#endif
};

void set_specific_boundary_values(Tensor& U, Tensor& V, int imax, int jmax) {
    //#ifdef _OMP
    // Outflow
    //#endif
#ifdef _OMP
#pragma omp parallel for
#endif
    for (int j=1; j!=jmax+1; ++j) {
        U(imax, j) = U(imax-1, j);
        V(imax+1, j) = V(imax, j);
    }
#ifdef _OMP
#pragma omp parallel for
#endif
    for (int i=1; i!=imax+1; ++i) {
        U(i, jmax+1) = U(i, jmax);
        V(i, jmax) = V(i, jmax-1);

        V(i, 0) = V(i, 1);
        U(i, 0) = U(i, 1);
    }
#ifdef _OMP
#pragma omp barrier
#endif

    // Left: input flow
    const float u = 0.8;
    const int width = 5;
#ifdef _OMP
#pragma omp parallel for
#endif
    for (int j=(jmax+1)/2-width/2; j!=(jmax+1)/2+width/2+1;++j)
    {
        U(0, j) = u;
    }
#ifdef _OMP
#pragma omp barrier
#endif
};

void compute_FG(Tensor& F, Tensor& G, const Tensor& U, const Tensor& V, const Domain& domain, float dt, float Re, float dx, float dy, float gamma, int imax, int jmax, float gx, float gy) {
    // F = u + dt * (1./Re * (dudxdx + dudydy) - duudx - duvdy + gx);       // i=1..imax-1 j=1..jmax
    // G = v + dt * (1./Re * (dvdxdx + dvdydy) - duvdx - dvvdy + gy);       // i=1..imax   j=1..jmax-1
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


#ifdef _OMP
#pragma omp parallel for default(none)          \
private(dudxdx, dudydy, duudx, duvdy, dvdxdx, dvdydy, duvdx, dvvdy, uijm, uij, uijp, uimj, uimjp, uipj, vijm, vij, vijp, vimj, vipjm, vipj, one, two, three, four, five, six)                         \
shared(F, G, U, V, imax, jmax, domain, dxinv2, dyinv2, dxinv4, dyinv4, gammadxinv4, gammadyinv4, dt, Reinv, gx, gy)
#endif
    for (int i=0; i!=imax+2; ++i) {
        for (int j=0; j!=jmax+2; ++j) {
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
                    dudxdx = (uipj - 2*uij + uimj) * dxinv2;
                    dudydy = (uijp - 2*uij + uijm) * dyinv2;

                    one   = uij  + uipj;
                    two   = uimj + uij;
                    three = uij  - uipj;
                    four  = uimj - uij;
                    duudx =        dxinv4 * (      one  * one   -      two  * two )       \
                            + gammadxinv4 * ( std::fabs(one) * three - std::fabs(two) * four );

                    one   = vij + vipj;
                    two   = uij + uijp;
                    three = vijm + vipjm;
                    four  = uijm + uij;
                    five = uij - uijp;
                    six =  uijm - uij;
                    duvdy =        dyinv4 * (      one  * two -        three  * four )      \
                            + gammadyinv4 * ( std::fabs(one) * five -  std::fabs(three) * six );

                    F(i, j) = uij + dt * (Reinv * (dudxdx + dudydy) - duudx - duvdy + gx);
                }

                if (j<jmax+1 & domain(i, j+1).obstacle == false)
                {
                    dvdxdx = (vipj - 2*vij + vimj) * dxinv2;
                    dvdydy = (vijp - 2*vij + vijm) * dyinv2;

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
#ifdef _OMP
#pragma omp barrier
#endif
};

void compute_rhs_pressure(Tensor& RHS, const Tensor& F, const Tensor& G, const Domain& domain, float dx, float dy, float dt, int imax, int jmax) {
    const static float dxinv = 1.0f/dx;
    const static float dyinv = 1.0f/dy;

    const float dtinv = 1.0f / dt;

    float fij;
    float fimj;
    float gij;
    float gijm;

#ifdef _OMP
#pragma omp parallel for default(none) private(fij, fimj, gij, gijm) shared(RHS, F, G, domain, imax, jmax, dxinv, dyinv, dtinv)
#endif
    for (int i=0; i<imax+2; ++i) {
        for (int j=0; j<jmax+2; ++j) {
            if (domain(i, j).obstacle==false) {
                fimj = F(i-1, j);
                fij  = F(i, j);
                gijm = G(i, j-1);
                gij  = G(i, j);
                RHS(i, j) = dtinv * ( dxinv * ( fij - fimj ) + dyinv * ( gij - gijm ) );
            }
        }
    }
#ifdef _OMP
#pragma omp barrier
#endif
};

void SOR(Tensor& P, const Tensor& RHS, const Domain& domain, float& rit, float omega, float dx, float dy, int imax, int jmax) {
    const static float dxinv = 1.0f / dx;
    const static float dyinv = 1.0f / dy;
    const static float dxinv2 = dxinv * dxinv;
    const static float dyinv2 = dyinv * dyinv;
    const static float coeff = omega / ( 2 * (dxinv2 + dyinv2) );


    float rit_tmp;

    rit = 0.0f;

    float tmp_p;
    int edges;

    // Set pressure boundary conditions
#ifdef _OMP
#pragma omp parallel for private(tmp_p, edges)
#endif
    for (int i=0; i!=imax+2; ++i) {
        for (int j=0; j!=jmax+2; ++j) {
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
    for (int i=1; i!=imax+1; ++i) {
        for (int j=1; j!=jmax+1; ++j) {
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
#ifdef _OMP
#pragma omp parallel for reduction(max: rit) shared(imax, jmax, P) private(rit_tmp)
#endif
    for (int i=1; i!=imax+1; ++i) {
        for (int j=1; j!=jmax+1; ++j) {
            if (domain(i, j).obstacle == false) {
                rit_tmp = dxinv2 * ( P(i-1, j) - 2.0f * P(i, j) + P(i+1, j) ) \
                        + dyinv2 * ( P(i, j-1) - 2.0f * P(i, j) + P(i, j+1) ) \
                        - RHS(i, j);
                rit = std::max(std::fabs(rit_tmp), rit);
            }
        }
    }
};

void compute_uv(Tensor& U, Tensor& V, const Tensor& F, const Tensor& G, const Tensor& P, const Domain& domain, float dx, float dy, float dt, int imax, int jmax) {
    const static float dxinv = 1.0f / dx;
    const static float dyinv = 1.0f / dy;

    float dtdx = dt * dxinv;
    float dtdy = dt * dyinv;

    float pipj;
    float pij;
    float pijp;

#ifdef _OMP
#pragma omp parallel for private(pij, pijp, pipj)
#endif
    for (int i=1; i!=imax+1; ++i) {
        for (int j=1; j!=jmax+1; ++j) {
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
