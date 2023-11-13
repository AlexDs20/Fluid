#pragma once
#include <chrono>
#include <iostream>

class Timer {
    public:
        Timer(int it=1) {
            iteration = it;
            start = std::chrono::system_clock::now();
        };
        ~Timer(){
            std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
            std::cout << "Elapsed time: "
                << (std::chrono::duration_cast<std::chrono::microseconds>(end - start)).count()/(float) iteration
                << "µs"
                << std::endl;
        };
    private:
        int iteration;
        std::chrono::time_point<std::chrono::system_clock> start;
};
