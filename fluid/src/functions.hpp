#pragma once
#include "tensor.hpp"

float adaptive_time_step_size( const Tensor& U, const Tensor& V, float dx, float dy, float Re, float tau, float dt=0.5f);

void set_boundary_values(Tensor& U, Tensor& V);

void set_object_boundary_values(Tensor& U, Tensor& V, Tensor& OBS);

void compute_FG(Tensor& F, Tensor& G, const Tensor& U, const Tensor& V, float dt, float Re, float dx, float dy, float gamma);

void compute_rhs_pressure(Tensor& RHS, const Tensor& F, const Tensor& G, float dx, float dy, float dt);

void SOR(Tensor& P, float& rit, const Tensor& RHS, float omega, float dx, float dy);

void compute_rit(float& rit);

void compute_uv();


