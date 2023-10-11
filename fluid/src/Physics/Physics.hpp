#pragma once
#include "parameters.hpp"
#include "Math/tensor.hpp"
#include "Message/Message.hpp"
#include "Physics/calculate.hpp"


class Fluid{
public:
    Fluid(){
        // Initialize the fields
        U = new Tensor({p.imax+2, p.jmax+2}, p.u0);
        V = new Tensor({p.imax+2, p.jmax+2}, p.v0);
        P = new Tensor({p.imax+2, p.jmax+2}, p.p0);
        F = new Tensor({p.imax+2, p.jmax+2}, 0.0f);
        G = new Tensor({p.imax+2, p.jmax+2}, 0.0f);
        RHS = new Tensor({p.imax+2, p.jmax+2}, 0.0f);
        domain = new Domain({p.imax+2, p.jmax+2}, Cell());

        // Set flags for edges and obstacle     (because obstacles don't move)
        set_constant_flags(*domain, p.imax, p.jmax);
    };
    ~Fluid(){
        delete U;
        delete V;
        delete P;
        delete F;
        delete G;
        delete RHS;
        delete domain;
    };

    // Run simulation until simulated long enough
    void update(){
        set_boundary_values(*U, *V, *domain, p.imax, p.jmax);
        set_specific_boundary_values(*U, *V, p.imax, p.jmax);
        float dt = adaptive_time_step_size(*U, *V, p.dx, p.dy, p.Re, p.tau, p.dt_max, p.imax, p.jmax);
        compute_FG(*F, *G, *U, *V, *domain, dt, p.Re, p.dx, p.dy, p.gamma, p.imax, p.jmax, p.gx, p.gy);
        compute_rhs_pressure(*RHS, *F, *G, *domain, p.dx, p.dy, dt, p.imax, p.jmax);

        int it = 0;
        float rit = 0;
        do {
            ++it;
            SOR(*P, *RHS, *domain, rit, p.omega, p.dx, p.dy, p.imax, p.jmax);
        } while (it < p.it_max && rit > p.eps);

        compute_uv(*U, *V, *F, *G, *P, *domain, p.dx, p.dy, dt, p.imax, p.jmax);
    };

public:
    Parameters p;
    Tensor *U;
    Tensor *V;
    Tensor *P;
    Tensor *F;
    Tensor *G;
    Tensor *RHS;
    Domain *domain;
};


class Physics: public Sender{
public:
    Physics(MessageBus* m): Sender(m) {};
    void initialise(){
    };
    void reset(){};
    void update(float dt){};
    static void read_message(Message){};
private:
};
