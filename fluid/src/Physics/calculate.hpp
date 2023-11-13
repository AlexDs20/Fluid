#pragma once
#include <string>
#include "Math/Matrix.hpp"

struct Constants {
    int imax;
    int jmax;
    float dx;
    float dy;
    float dxinv;
    float dyinv;
    float dxdxinv;
    float dydyinv;
    float dxinv4;
    float dyinv4;
    float gammadxinv4;
    float gammadyinv4;
    float Re;
    float tau;
    float omega;
    float coeff;
    // float Re_dt;         // need to remove 1 to fit 1 cache line
};

struct Parameters {
    // Fluid
    float Re;
    float u0;
    float v0;
    float p0;
    float gx;
    float gy;

    // time
    float t_max;
    float dt_max;
    float tau;

    // grid
    int imax;
    int jmax;
    float xlength;
    float ylength;

    // pressure
    int it_max;
    float eps;
    float omega;

    // Derivative scheme
    float gamma;
};

void set_constant_flags(Matrixi&, int imax, int jmax);

void set_boundary_values(Matrix& U, Matrix& V, const Matrixi&, int imax, int jmax);

void set_specific_boundary_values(Matrix& U, Matrix& V, int imax, int jmax);

float adaptive_time_step_size(Matrix& U, Matrix& V, float dt, const Constants&);

void compute_FG(Matrix& F, Matrix& G, const Matrix& U, const Matrix& V, const Matrixi& domain, float dt, float gx, float gy, const Constants& cst);

void compute_rhs_pressure(Matrix& RHS, Matrix& F, Matrix& G, Matrixi& domain, float dt, const Constants& c);

void SOR(Matrix& P, const Matrix& RHS, const Matrixi&, float& rit, const Constants& p);

void compute_uv(Matrix& U, Matrix& V, Matrix& F, Matrix& G, Matrix& P, Matrixi&, float dt, const Constants& cst);

void get_parameters(const std::string& problem, Parameters& params, Constants& constants);
