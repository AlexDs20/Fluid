#pragma once
#include "Message/Message.hpp"
#include "Renderer/object.hpp"

class Scene: public Sender {
public:
    Scene(MessageBus* mBus): Sender(mBus) {
    };

    void read_message(Message m){
    };
private:
};
