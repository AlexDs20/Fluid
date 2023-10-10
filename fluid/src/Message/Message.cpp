/*
#include <iostream>
#include <vector>
#include <queue>
#include "MessageBus.hpp"


class Phy : public Sender {
    public:
        static void read_message(Message m){
            std::cout << "Reading message from phy: " << m.e << std::endl;
        };
        Phy(MessageBus* bus): Sender(bus){
        };
};


int main(){
    MessageBus B;
    Phy phy(&B);
    B.add_receiver(phy.read_message);

    Message m{Event::KEY_PRESS};
    phy.send_message(m);
    B.dispatch();
};
*/
