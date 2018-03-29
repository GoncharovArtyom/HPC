//
// Created by artyom on 26.03.18.
//

#ifndef AUTOREG_TIMER_H
#define AUTOREG_TIMER_H

#include <fstream>
#include <string>
#include <ctime>

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
            start = clock();
        }

        void end_clock() {
            end = clock();
        }

        void write(std::string function_name) {
            std::clock_t diff = end - start;
            double n_secs = double(diff) / CLOCKS_PER_SEC;

            out <<first_second_columns<< function_name << "\t\t" << n_secs <<" sec"<< std::endl;
        }

    private:
        std::ofstream out;
        std::clock_t start;
        std::clock_t end;
        std::string first_second_columns;
    };
}
#endif //AUTOREG_TIMER_H
