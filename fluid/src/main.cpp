#include <iostream>
#include <string>
#include <math.h>

#include "tensor.hpp"
#include "functions.hpp"
#include "parameters.hpp"

int main(int argc, char** argv) {
    std::string config;
    if (argc>1)
        config = argv[1];

    Parameters p;

    // Initialize the fields
    // TODO: Investigate if order of the dim / grid or if grouping the fields together changes the perf.
    Tensor U({p.imax+2, p.jmax+2}, p.u0);
    Tensor V({p.imax+2, p.jmax+2}, p.v0);
    Tensor P({p.imax+2, p.jmax+2}, p.p0);
    Tensor EXT_X({p.imax+2, p.jmax+2});     // External forces
    Tensor EXT_Y({p.imax+2, p.jmax+2});     // External forces
    Tensor F({p.imax+2, p.jmax+2});
    Tensor G({p.imax+2, p.jmax+2});
    Tensor RHS({p.imax+2, p.jmax+2});
    Tensor OBS({p.imax+2, p.jmax+2});

    while (p.t < p.t_max) {
        std::cout << p.t << "/" << p.t_max << "\t" << p.dt << std::endl;
        p.dt = adaptive_time_step_size(U, V, p.dx, p.dy, p.Re, p.tau, p.dt);
        set_boundary_values(U, V);
        set_object_boundary_values(U, V, OBS);
        compute_FG(F, G, U, V, p.dt, p.Re, p.dx, p.dy, p.gamma);
        compute_rhs_pressure(RHS, F, G, p.dx, p.dy, p.dt);

        p.it = 0;
        while (p.it<p.it_max && p.rit > p.eps * p.norm_p0) {
            SOR(P, p.rit, RHS, p.omega, p.dx, p.dy);
            ++p.it;
        }

        compute_uv(U, V, F, G, P, p.dx, p.dy, p.dt);

        p.t += p.dt;
        p.n += 1;
    }

    return 0;
}
