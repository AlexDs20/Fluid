#pragma once
#include "Math/tensor.hpp"

struct Boundary {
    Tensor cellType;
    Tensor NCell;
    Tensor SCell;
    Tensor WCell;
    Tensor ECell;
    Tensor boundaryType;
    Tensor test;
};

float adaptive_time_step_size( const Tensor& U, const Tensor& V, float dx, float dy, float Re, float tau, float dt, int imax, int jmax);

void set_constant_flags(Boundary& boundary, int imax, int jmax);

void set_boundary_values(Tensor& U, Tensor& V, const Boundary& boundary, int imax, int jmax);

void set_specific_boundary_values(Tensor& U, Tensor& V, int imax, int jmax);

void compute_FG(Tensor& F, Tensor& G, const Tensor& U, const Tensor& V, const Boundary&, float dt, float Re, float dx, float dy, float gamma, int imax, int jmax);

void compute_rhs_pressure(Tensor& RHS, const Tensor& F, const Tensor& G, const Boundary&, float dx, float dy, float dt, int imax, int jmax);

void SOR(Tensor& P, const Tensor& RHS, const Boundary&, float& rit, float omega, float dx, float dy, int imax, int jmax);

void compute_uv(Tensor& U, Tensor& V, const Tensor& F, const Tensor& G, const Tensor& P, const Boundary&, float dx, float dy, float dt, int imax, int jmax);
