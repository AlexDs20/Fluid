#include <iostream>
#include <cmath>
#include "Physics/calculate.hpp"


float adaptive_time_step_size( const Tensor& U, const Tensor& V, float dx, float dy, float Re, float tau, float dt) {
    const static float dxinv2 = 1.0f / (dx*dx);
    const static float dyinv2 = 1.0f / (dy*dy);

    const static float Re_dt = 0.5f * Re / (dxinv2 + dyinv2);

    if ( (tau <= 0) | (tau > 1) ) {
        return dt;
    }

    float absumax=0, absvmax=0;
    for (int i=0; i<U.imax()+2; ++i){
        for (int j=0; j<U.jmax()+2; ++j) {
            absumax = std::max(absumax, std::fabs(U({i, j})));
            absvmax = std::max(absvmax, std::fabs(V({i, j})));
        }
    }

    const float dt_v_max = std::min(dx / absumax, dy / absvmax);
    // const float dt_v_max = std::min(dx / U.amax(), dy / V.amax());

    float ret = std::min(Re_dt, dt_v_max);
    ret = std::min(dt, tau*ret);

    return ret;
};


void set_constant_flags(Boundary& boundary) {
    // Set the flags that do not change over time
    // Boundaries and sphere in the middle
    const static int imax = boundary.cellType.imax();
    const static int jmax = boundary.cellType.jmax();
    const static float isBoundary = 1;
    const static float isFluid = 0;

    // Sphere info
    const int px = (imax+2) / 3;
    const int py = (jmax+2) / 2;
    const float radius = std::min(imax, jmax) / 10;

    for (int i=0; i<imax+2; ++i) {
        for (int j=0; j<jmax+2; ++j) {
            // Above and Below
            if (i==0 || i==imax+1) {
                boundary.cellType({i, j}) = isBoundary;
            } else if (j==0 || j==jmax+1) {      // left / right
                boundary.cellType({i, j}) = isBoundary;
            } else if ( (px-i)*(px-i) + (py-j)*(py-j) < radius*radius) {    // close to sphere
                boundary.cellType({i, j}) = isBoundary;
            } else {
                boundary.cellType({i, j}) = isFluid;
            }
        }
    }

    for (int i=0; i<imax+2; ++i) {
        for (int j=0; j<jmax+2; ++j) {
            if (boundary.cellType({i, j}) == isBoundary) {
                if (i<imax+1 & boundary.cellType({i+1, j}) == isFluid){
                    boundary.ECell({i, j}) = 1;
                }
                if (i>0 & boundary.cellType({i-1, j}) == isFluid){
                    boundary.WCell({i, j}) = 1;
                }
                if (j<jmax+1 & boundary.cellType({i, j+1}) == isFluid){
                    boundary.NCell({i, j}) = 1;
                }
                if (j>0 & boundary.cellType({i, j-1}) == isFluid){
                    boundary.SCell({i, j}) = 1;
                }
            }
        }
    }
};


void set_boundary_flags(Boundary& boundary) {
    const static int imax = boundary.cellType.imax();
    const static int jmax = boundary.cellType.jmax();
    const static float isBoundary = 1.0f;
    const static float isFluid = 0.0f;

    // Put Sphere
    int px = (imax+2) / 3;
    int py = (jmax+2) / 2;
    float radius = std::min(imax, jmax) / 10;

    // Set all grid cells of obstacles as boundary
    for (int i=0; i<imax+2; ++i) {
        // boundary.cellType({i, 0}) = isBoundary;
        // boundary.cellType({i, jmax+1}) = isBoundary;

        for (int j=0; j<jmax+2; ++j) {
            if (i==0 ||  i==(imax+1) || j==0 || j==(jmax+1))
                boundary.cellType({i, j}) = isBoundary;
            else
                boundary.cellType({i, j}) = isFluid;

            // boundary.cellType({0, j}) = isBoundary;
            // boundary.cellType({imax+1, j}) = isBoundary;

            // if ( (px-i)*(px-i)+(py-j)*(py-j) < radius*radius )
            //     boundary.cellType({i, j}) = isBoundary;
        }
    }

    for (int i=0; i<imax+2; ++i) {
        for (int j=0; j<jmax+2; ++j) {
            if (boundary.cellType({i, j})==isBoundary) {
                if ((j!=jmax+1) && (boundary.cellType({i, j+1}) == isFluid))
                    boundary.NCell({i, j}) = 1.0f;
                if ((j!=0) && (boundary.cellType({i, j-1}) == isFluid))
                    boundary.SCell({i, j}) = 1.0f;
                if ((i!=0) && (boundary.cellType({i-1, j}) == isFluid))
                    boundary.WCell({i, j}) = 1.0f;
                if ((i!=imax+1) && (boundary.cellType({i+1, j}) == isFluid))
                    boundary.ECell({i, j}) = 1.0f;
            }
        }
    }

    // TODO: Check if we have boundaries with fluid on opposite sides or 3 sides
};


void set_boundary_values(Tensor& U, Tensor& V, const Boundary& boundary){
    // Currently hard-code up/down left/right
    const static int imax = U.imax();
    const static int jmax = U.jmax();
    const static float isBoundary = 1;
    const static float isFluid = 0;


    for (int i=0; i<imax+2; ++i) {
        for (int j=0; j<jmax+2; ++j) {
            if (boundary.cellType({i, j}) == isBoundary) {
                // TODO: check boundary type
                //       Start with no slip
                // U
                if (boundary.WCell({i, j})) {
                    U({i-1, j}) = 0;
                } else {
                    if (boundary.NCell({i, j})) {
                        U({i-1, j}) = -U({i-1, j+1});
                    }
                    if (boundary.SCell({i, j})) {
                        U({i-1, j}) = -U({i-1, j-1});
                    }
                }
                if (boundary.ECell({i, j})) {
                    U({i, j}) = 0;
                } else {
                    if (boundary.NCell({i, j})) {
                        U({i, j}) = -U({i, j+1});
                    }
                    if (boundary.SCell({i, j})) {
                        U({i, j}) = -U({i, j-1});
                    }
                }

                // V
                if (boundary.NCell({i, j})) {
                    V({i, j}) = 0;
                } else {
                    if (boundary.WCell({i, j})) {
                        V({i, j}) = -V({i-1, j});
                    }
                    if (boundary.ECell({i, j})) {
                        V({i, j}) = -V({i+1, j});
                    }
                }
                if (boundary.SCell({i, j})) {
                    V({i, j-1}) = 0;
                } else {
                    if (boundary.WCell({i, j})) {
                        V({i, j-1}) = -V({i-1, j-1});
                    }
                    if (boundary.ECell({i, j})) {
                        V({i, j-1}) = -V({i+1, j-1});
                    }
                }
            }
        }
    }
};

void set_specific_boundary_values(Tensor& U, Tensor& V) {
    const static int imax = U.imax();
    const static int jmax = U.jmax();

    // Right: outflow
    for (int j=1; j!=jmax+1; ++j) {
        U({imax, j}) = U({imax-1, j});
        V({imax+1, j}) = V({imax, j});
    }
    // for (int i=1; i!=imax+1; ++i) {
    //     U({i, jmax+1}) = U({i, jmax});
    //     V({i, jmax}) = V({i, jmax-1});
    // }

    // Left: input flow
    const float u = 0.8;
    const int width = 5;
    for (int j=(jmax+1)/2-width/2; j!=(jmax+1)/2+width/2+1;++j)
    {
        U({0, j}) = u;
    }

    // const float v = 0.8;
    // for (int i=(imax+1)/2-width/2; i!=(imax+1)/2+width/2+1;++i)
    // {
    //     V({i, 0}) = v;
    // }
};

void compute_FG(Tensor& F, Tensor& G, const Tensor& U, const Tensor& V, float dt, float Re, float dx, float dy, float gamma, const Boundary& boundary) {
    // F = u + dt * (1./Re * (dudxdx + dudydy) - duudx - duvdy + gx);       // i=1..imax-1 j=1..jmax
    // G = v + dt * (1./Re * (dvdxdx + dvdydy) - duvdx - dvvdy + gy);       // i=1..imax   j=1..jmax-1

    const static float isBoundary = 1.0f;
    const static float isFluid = 0.0f;

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

    // TODO
    const static float gx = 0.0f, gy = 0.0f;

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


    for (int i=0; i!=imax+2; ++i) {
        for (int j=0; j!=jmax+2; ++j) {
            if (boundary.cellType({i, j}) == isFluid) {
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

                if (i<imax+1 & boundary.cellType({i+1, j}) == isFluid)
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

                    F({i, j}) = uij + dt * (Reinv * (dudxdx + dudydy) - duudx - duvdy + gx);
                }

                if (j<jmax+1 & boundary.cellType({i, j+1})==isFluid)
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

                    G({i, j}) = vij + dt * (Reinv * (dvdxdx + dvdydy) - duvdx - dvvdy + gy );

                }
            } else {
                if (boundary.NCell({i, j}))
                    G({i, j}) = V({i, j});
                else if (boundary.SCell({i, j}))
                    G({i, j-1}) = V({i, j-1});

                if (boundary.WCell({i, j}))
                    F({i-1, j}) = U({i-1, j});
                else if (boundary.ECell({i, j}))
                    F({i, j}) = U({i, j});
            }
        }
    }
};

void compute_rhs_pressure(Tensor& RHS, const Tensor& F, const Tensor& G, float dx, float dy, float dt, const Boundary& boundary) {
    const static float dxinv = 1.0f/dx;
    const static float dyinv = 1.0f/dy;
    const static float imax = F.imax();
    const static float jmax = F.jmax();

    const float dtinv = 1.0f / dt;

    float fij;
    float fimj;
    float gij;
    float gijm;

    const static float isBoundary = 1.0f;
    const static float isFluid = 0.0f;

    for (int i=0; i<imax+2; ++i) {
        for (int j=0; j<jmax+2; ++j) {
            if (boundary.cellType({i, j})==isFluid) {       // TODO: Could delete this?
                fimj = F({i-1, j});
                fij  = F({i, j});
                gijm = G({i, j-1});
                gij  = G({i, j});
                RHS({i, j}) = dtinv * ( dxinv * ( fij - fimj ) + dyinv * ( gij - gijm ) );
            }
        }
    }
};

void SOR(Tensor& P, float& rit, const Tensor& RHS, float omega, float dx, float dy, const Boundary& boundary) {
    const static int imax = P.imax();
    const static int jmax = P.jmax();
    const static float dxinv = 1.0f / dx;
    const static float dyinv = 1.0f / dy;
    const static float dxinv2 = dxinv * dxinv;
    const static float dyinv2 = dyinv * dyinv;
    const static float coeff = omega / ( 2 * (dxinv2 + dyinv2) );


    float rit_tmp;

    rit = 0.0f;

    float tmp_p;
    int edges;
    const static float isBoundary = 1.0f;
    const static float isFluid = 0.0f;

    // Set pressure boundary conditions
    for (int i=0; i!=imax+2; ++i) {
        for (int j=0; j!=jmax+2; ++j) {
            if (boundary.cellType({i, j})==isBoundary) {
                tmp_p = 0.0f;
                edges = 0;
                if (boundary.NCell({i, j})) {
                    tmp_p += P({i, j+1});
                    ++edges;
                }
                if (boundary.SCell({i, j})) {
                    tmp_p += P({i, j-1});
                    ++edges;
                }
                if (boundary.ECell({i, j})) {
                    tmp_p += P({i+1, j});
                    ++edges;
                }
                if (boundary.WCell({i, j})) {
                    tmp_p += P({i-1, j});
                    ++edges;
                }
                P({i, j}) = tmp_p / (float)edges;
            }
        }
    }


    for (int i=1; i!=imax+1; ++i) {
        for (int j=1; j!=jmax+1; ++j) {
            if (boundary.cellType({i, j}) == isFluid) {
                P({i, j}) = (1.0f - omega) * P({i, j})                  \
                    + coeff                                             \
                    * ( dxinv2 * ( P({i+1, j}) + P({i-1, j}) )          \
                      + dyinv2 * ( P({i, j+1}) + P({i, j-1}) )          \
                      - RHS({i, j})                                     \
                      );
            }
        }
    }

    for (int i=1; i!=imax+1; ++i) {
        for (int j=1; j!=jmax+1; ++j) {
            if (boundary.cellType({i, j}) == isFluid) {
                rit_tmp = dxinv2 * ( P({i-1, j}) - 2.0f * P({i, j}) + P({i+1, j}) ) \
                        + dyinv2 * ( P({i, j-1}) - 2.0f * P({i, j}) + P({i, j+1}) ) \
                        - RHS({i, j});
                rit = std::max(std::fabs(rit_tmp), rit);
            }
        }
    }
};

void compute_uv(Tensor& U, Tensor& V, const Tensor& F, const Tensor& G, const Tensor& P, float dx, float dy, float dt, const Boundary& boundary) {
    const static float dxinv = 1.0f / dx;
    const static float dyinv = 1.0f / dy;

    const static float imax = U.imax();
    const static float jmax = U.jmax();

    const static float isBoundary = 1.0f;
    const static float isFluid = 0.0f;

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

            if ( (boundary.cellType({i, j})==isFluid) && boundary.cellType({i+1, j})==isFluid )
                U({i, j}) = F({i, j}) - dtdx * (pipj - pij);
            if ( (boundary.cellType({i, j})==isFluid) && boundary.cellType({i, j+1})==isFluid )
                V({i, j}) = G({i, j}) - dtdy * (pijp - pij);
        }
    }
};
