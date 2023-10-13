#pragma once
#include "Console/Console.hpp"
#include "Message/Message.hpp"
#include <ctime>
#include <iostream>
#include <chrono>


class Console{
public:
    static void read_message(Message m){
        const auto now = std::chrono::system_clock::now();
        const std::time_t t_c = std::chrono::system_clock::to_time_t(now);

        char buffer[100];
        if (std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", std::localtime(&t_c)))
            std::cout << "Console: " << buffer << ": " << m.e << std::endl;
    };
};
