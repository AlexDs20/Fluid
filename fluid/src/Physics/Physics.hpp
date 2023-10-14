#pragma once
#include "parameters.hpp"
#include "Math/tensor.hpp"
#include "Message/Message.hpp"
#include "Physics/calculate.hpp"

#include "utils.hpp"

#include <cstdlib>
#include <time.h>
using namespace std;

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
    // {
    //     Timer t(1000);
    // }
    void update(){
        set_boundary_values(*U, *V, *domain, p.imax, p.jmax);                   // 66 -> 120 mus

        set_specific_boundary_values(*U, *V, p.imax, p.jmax);                   // 4 -> 7 mus
        float dt = adaptive_time_step_size(*U, *V, p.dt_max, p);                // 95 -> 125 mus
        compute_FG(*F, *G, *U, *V, *domain, dt, p);                             // 1100 -> 1500
        compute_rhs_pressure(*RHS, *F, *G, *domain, dt, p);                     // 240 -> 320 mus

        int it_max = 20;
        int it = 0;
        float rit = 0;
        do {
            ++it;
            SOR(*P, *RHS, *domain, rit, p);                                     // 840 -> 1150 mus
        } while (it < it_max && rit > p.eps);

        compute_uv(*U, *V, *F, *G, *P, *domain, dt, p);                         // 400 -> 670 mus
    };

public:
    __attribute__ ((aligned(64))) Parameters p;
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
