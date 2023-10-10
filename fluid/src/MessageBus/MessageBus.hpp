#pragma once
#include <vector>
#include <queue>

enum Event {
    KEY_PRESS,
};


struct Message{
        Event e;
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
