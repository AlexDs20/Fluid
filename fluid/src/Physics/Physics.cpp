#include <string>
#include "Math/Matrix.hpp"
#include "Message/Message.hpp"
#include "Physics/calculate.hpp"
#include "Physics/Physics.hpp"
#include "utils.hpp"


Fluid::Fluid(const std::string& problem, MessageBus* m): Sender(m) {
    get_parameters(problem, params, constants);
    // Initialize the fields
    U = new Matrix(params.imax+2, params.jmax+2, params.u0);
    V = new Matrix(params.imax+2, params.jmax+2, params.v0);
    P = new Matrix(params.imax+2, params.jmax+2, params.p0);
    F = new Matrix(params.imax+2, params.jmax+2, 0.0f);
    G = new Matrix(params.imax+2, params.jmax+2, 0.0f);
    RHS = new Matrix(params.imax+2, params.jmax+2, 0.0f);
    domain = new Matrixi(params.imax+2, params.jmax+2, 0);

    // Set flags for edges and obstacle     (because obstacles don't move)
    set_constant_flags(*domain, params.imax, params.jmax);                          // probably need the std::string problem here
                                                                                    // Scalar: 220 -> 430 mus
                                                                                    // Vector: 100 mus
};

Fluid::~Fluid(){
    delete U;
    delete V;
    delete P;
    delete F;
    delete G;
    delete RHS;
    delete domain;
};

// Run simulation until simulated long enough
// { int n = 200; Timer t(n); for (int it=0; it<n; ++it)
// }
void Fluid::update(){
    float dt;
    set_boundary_values(*U, *V, *domain, constants.imax, constants.jmax);            // 70 -- 81 mus
    set_specific_boundary_values(*U, *V, constants.imax, constants.jmax);            // 4 mus

    dt = adaptive_time_step_size(*U, *V, params.dt_max, constants);                  // 95 -- 125 mus
                                                                                     // Vector: 25 -> 50
    compute_FG(*F, *G, *U, *V, *domain, dt, params.gx, params.gy, constants);        // 1100 -- 1500
    compute_rhs_pressure(*RHS, *F, *G, *domain, dt, constants);                      // 182 -- 215 mus
                                                                                     // Vector: 75 -- 150 mus

    int it_max = 5;
    int it = 0;
    float rit = 0;
    do {
        ++it;
        SOR(*P, *RHS, *domain, rit, constants);                                      // 840 -- 1150 mus
    } while (it < it_max && rit > params.eps);

    compute_uv(*U, *V, *F, *G, *P, *domain, dt, constants);                          // 400 -- 670 mus
                                                                                     // Vector: 125 -- 310
};
