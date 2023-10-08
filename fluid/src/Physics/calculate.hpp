#pragma once
#include "Math/tensor.hpp"


float adaptive_time_step_size( const Tensor& U, const Tensor& V, float dx, float dy, float Re, float tau, float dt, int imax, int jmax);

void set_constant_flags(Domain&, int imax, int jmax);

void set_boundary_values(Tensor& U, Tensor& V, const Domain&, int imax, int jmax);

void set_specific_boundary_values(Tensor& U, Tensor& V, int imax, int jmax);

void compute_FG(Tensor& F, Tensor& G, const Tensor& U, const Tensor& V, const Domain&, float dt, float Re, float dx, float dy, float gamma, int imax, int jmax);

void compute_rhs_pressure(Tensor& RHS, const Tensor& F, const Tensor& G, const Domain&, float dx, float dy, float dt, int imax, int jmax);

void SOR(Tensor& P, const Tensor& RHS, const Domain&, float& rit, float omega, float dx, float dy, int imax, int jmax);

void compute_uv(Tensor& U, Tensor& V, const Tensor& F, const Tensor& G, const Tensor& P, const Domain&, float dx, float dy, float dt, int imax, int jmax);
