#pragma once
#include "Math/tensor.hpp"
#include "parameters.hpp"

float adaptive_time_step_size( const Tensor& U, const Tensor& V, float dt, const Parameters&);

void set_constant_flags(Domain&, int imax, int jmax);

void set_boundary_values(Tensor& U, Tensor& V, const Domain&, int imax, int jmax);

void set_specific_boundary_values(Tensor& U, Tensor& V, int imax, int jmax);

void compute_FG(Tensor& F, Tensor& G, const Tensor& U, const Tensor& V, const Domain& domain, float dt, const Parameters& p);

void compute_rhs_pressure(Tensor& RHS, const Tensor& F, const Tensor& G, const Domain& domain, float dt, const Parameters& p);

void SOR(Tensor& P, const Tensor& RHS, const Domain&, float& rit, const Parameters& p);

void compute_uv(Tensor& U, Tensor& V, const Tensor& F, const Tensor& G, const Tensor& P, const Domain&, float dt, const Parameters&);
