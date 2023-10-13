#pragma once
#include <GLFW/glfw3.h>
#include "Message/Message.hpp"

class Input: public Sender{
public:
    Input(MessageBus *mBus, GLFWwindow* window): Sender(mBus), window(window){
    };
    void update(){
        if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
            send_message(Message {Event::KEY_PRESS_ESCAPE});
        if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
            send_message(Message {Event::KEY_PRESS_W});
        if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
            send_message(Message {Event::KEY_PRESS_A});
        if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
            send_message(Message {Event::KEY_PRESS_S});
        if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
            send_message(Message {Event::KEY_PRESS_D});
    };

    static void read_message(Message m){
    };
private:
    GLFWwindow* window;
};
