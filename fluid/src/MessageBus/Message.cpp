#include <iostream>
#include <vector>
#include <queue>


class Message{
    public:
        int a = 1;
};

class MessageBus{
    public:
        void dispatch(){
            while(!messages.empty()) {
                for (auto& reader: readers){
                    reader(messages.front());
                }
                messages.pop();
            }
        };

        void add_receiver(void func(Message)){
            readers.push_back(func);
        };

        void send_message(Message m){
            messages.push(m);
        };
    private:
        // Need a list of receiver or receiver functions that can process a message
        std::vector<void (*)(Message)> readers;
        std::queue<Message> messages;
};

class Phy{
    public:
        static void read_message(Message m){
            std::cout << "Reading message from phy: " << m.a << std::endl;
        };
        Phy(){
        };
        Phy(MessageBus* bus): bus(bus){
        };
        void send_message(Message m){
            bus->send_message(m);
        };
    private:
        MessageBus* bus;
};


/*
int main(){
    MessageBus B;
    Phy phy(&B);
    B.add_receiver(phy.read_message);

    phy.send_message(Message());
    B.dispatch();
};
*/
