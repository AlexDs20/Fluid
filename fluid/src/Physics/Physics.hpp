#pragma once
#include <string>
#include "Math/tensor.hpp"
#include "Message/Message.hpp"
#include "Physics/calculate.hpp"


class Fluid: public Sender {
public:
    Fluid(const std::string& problem, MessageBus* m);
    ~Fluid();
    void update();

public:
    Parameters params;
    Constants constants;
    Tensor *U;
    Tensor *V;
    Tensor *P;
    Tensor *F;
    Tensor *G;
    Tensor *RHS;
    Domain *domain;
};

