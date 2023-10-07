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

void compute_FG(Tensor& F, Tensor& G, const Tensor& U, const Tensor& V, float dt, float Re, float dx, float dy, float gamma, const Boundary&, int imax, int jmax);

void compute_rhs_pressure(Tensor& RHS, const Tensor& F, const Tensor& G, float dx, float dy, float dt, const Boundary&, int imax, int jmax);

void SOR(Tensor& P, float& rit, const Tensor& RHS, float omega, float dx, float dy, const Boundary&, int imax, int jmax);

void compute_uv(Tensor& U, Tensor& V, const Tensor& F, const Tensor& G, const Tensor& P, float dx, float dy, float dt, const Boundary&, int imax, int jmax);
