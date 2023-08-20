#include <iostream>
#include <string>
#include <math.h>

#include "tensor.hpp"
#include "functions.hpp"


int main(int argc, char** argv) {
    std::string config;
    if (argc>1)
        config = argv[1];

    // Fluid
    float Re = 100.0f;
    float u0 = 3.0f;
    float v0 = 5.0f;
    float p0 = 0.0f;

    // time
    float t = 0;
    float t_max = 100;
    float dt=0.5;
    float tau = 0.5;

    // grid
    int dims = 2;
    int length_x = 2;
    int length_y = 1;
    int imax = 10;
    int jmax = 5;
    float dx = (float)length_x / imax;
    float dy = (float)length_y / jmax;

    // pressure
    int it = 0;
    int it_max = 10;
    float rit = 0.0f;
    float eps = 0.01f;
    float omega = 0.7;
    float gamma = 0.5;
    float norm_p0 = 0.5f;

    // Other stuff
    int n = 0;

    // Initialize the fields
    // TODO: Investigate if order of the dim / grid or if grouping the fields together changes the perf.
    Tensor U({imax+2, jmax+2}, u0);
    Tensor V({imax+2, jmax+2}, v0);
    Tensor P({imax+2, jmax+2}, p0);
    Tensor EXT_X({imax+2, jmax+2});     // External forces
    Tensor EXT_Y({imax+2, jmax+2});     // External forces
    Tensor F({imax+2, jmax+2});
    Tensor G({imax+2, jmax+2});
    Tensor RHS({imax+2, jmax+2});
    Tensor OBS({imax+2, jmax+2});

    while (n<3) {
        dt = adaptive_time_step_size(U, V, dx, dy, Re, tau, dt);
        set_boundary_values(U, V);
        set_object_boundary_values(U, V, OBS);
        compute_FG(F, G, U, V, dt, Re, dx, dy, gamma);
        compute_rhs_pressure(RHS, F, G, dx, dy, dt);

        std::cout << t << "/" << t_max << "\t" << dt << std::endl;

        it = 0;
        while (it<it_max && rit > eps * norm_p0) {
            SOR(P, rit, RHS, omega, dx, dy);
            //compute_rit(rit);
            ++it;
        }

        compute_uv();

        t += dt;
        n += 1;
    }

    return 0;
}
