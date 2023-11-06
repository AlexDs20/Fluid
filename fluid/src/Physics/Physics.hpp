#pragma once
#include <string>
#include "Math/Matrix.hpp"
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
    Matrix *U;
    Matrix *V;
    Matrix *P;
    Matrix *F;
    Matrix *G;
    Matrix *RHS;
    Matrixi *domain;
};

