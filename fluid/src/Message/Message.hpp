#pragma once
#include <vector>
#include <queue>

enum Event {
    KEY_PRESS_W,
    KEY_PRESS_A,
    KEY_PRESS_S,
    KEY_PRESS_D,
    KEY_PRESS_ESCAPE,

    SET_IMAX,
    SET_JMAX,
};

struct Message{
        Event e;
        int value;
};

class MessageBus{
    public:
        void dispatch(){
            while(!messages.empty()) {
                for (auto& receiver: receivers){
                    receiver(messages.front());
                }
                messages.pop();
            }
        };

        void add_receiver(void func(Message)){
            receivers.push_back(func);
        };

        void send_message(Message m){
            messages.push(m);
        };
    private:
        // Need a list of receiver functions that can process a message
        std::vector<void (*)(Message)> receivers;
        std::queue<Message> messages;
};

class Sender {
    public:
        Sender(MessageBus* m) : messageBus(m){};
        void send_message(Message m){
            messageBus->send_message(m);
        };
    protected:
        MessageBus* messageBus;
};
