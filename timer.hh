//
// Created by artyom on 26.03.18.
//

#ifndef AUTOREG_TIMER_H
#define AUTOREG_TIMER_H

#include <fstream>
#include <string>
#include <chrono>

namespace autoreg {
    struct FormattedTimer {
        FormattedTimer(std::string out_file_name, size3 const &zsize, size3 const &
        acf_size) :out(out_file_name) {
            std::string first_row = "zsize\t\tacf_size\tфункция\t\t\t\tвремя работы\n";
            out << first_row;

            first_second_columns = "(" + std::to_string(zsize(0)) + ", " + std::to_string(zsize(1))
                                   + ", " + std::to_string(zsize(2))
                                   + ")\t(" + std::to_string(acf_size(0)) + ", "
                                   + std::to_string(acf_size(1))
                                   + ", " + std::to_string(acf_size(2)) + ")\t";
        }

        void begin_clock() {
            start = std::chrono::system_clock::now();
        }

        void end_clock() {
            end = std::chrono::system_clock::now();
        }

        void write(std::string function_name) {
            auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
            out <<first_second_columns<< function_name << "\t\t" << diff <<" millisec"<< std::endl;
        }

    private:
        std::ofstream out;
        std::chrono::system_clock::time_point start;
        std::chrono::system_clock::time_point end;
        std::string first_second_columns;
    };
}
#endif //AUTOREG_TIMER_H
