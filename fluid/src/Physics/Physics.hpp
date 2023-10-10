#pragma once
#include "Message/Message.hpp"

class Physics: public Sender{
    public:
        Physics(MessageBus* m): Sender(m) {};
        void initialise(){};
        void reset(){};
        void update(float dt){};
        static void read_message(Message){};
};
