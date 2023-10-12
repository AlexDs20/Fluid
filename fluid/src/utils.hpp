#pragma once
#include <bits/chrono.h>
#include <chrono>
#include <iostream>

class Timer {
    public:
        Timer() {
            start = std::chrono::system_clock::now();
        };
        ~Timer(){
            std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
            std::cout << "Elapsed time: "
                << (std::chrono::duration_cast<std::chrono::microseconds>(end - start)).count()
                << "mus"
                << std::endl;
        };
    private:
        std::chrono::time_point<std::chrono::system_clock> start;
};
