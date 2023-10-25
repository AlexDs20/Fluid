#pragma once
#include <string>
#include "Math/tensor.hpp"

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

void set_constant_flags(Domain&, int imax, int jmax);

void set_boundary_values(Tensor& U, Tensor& V, const Domain&, int imax, int jmax);

void set_specific_boundary_values(Tensor& U, Tensor& V, int imax, int jmax);

float adaptive_time_step_size( const Tensor& U, const Tensor& V, float dt, const Constants&);

void compute_FG(Tensor& F, Tensor& G, const Tensor& U, const Tensor& V, const Domain& domain, float dt, float gx, float gy, const Constants& cst);

void compute_rhs_pressure(Tensor& RHS, const Tensor& F, const Tensor& G, const Domain& domain, float dt, const Constants& c);

void SOR(Tensor& P, const Tensor& RHS, const Domain&, float& rit, const Constants& p);

void compute_uv(Tensor& U, Tensor& V, const Tensor& F, const Tensor& G, const Tensor& P, const Domain&, float dt, const Constants& cst);

void get_parameters(const std::string& problem, Parameters& params, Constants& constants);
