#ifndef HPC_GENERATE_ZETA_PARALLEL_HH
#define HPC_GENERATE_ZETA_PARALLEL_HH

#include <vector>
#include <mutex>
#include <list>
#include <blitz/array.h>
#include <cmath>
#include <algorithm>

namespace parallel {
    using namespace std;

    struct ZetaGenerationBlock {
        int t_start;
        int x_start;
        int y_start;

        int t_end;
        int x_end;
        int y_end;

        int x_id;
        int y_id;
        int t_id;
    };

    struct ZetaGenerationController {

        ZetaGenerationController(int t_step_, int x_step_, int y_step_, int t_max, int x_max, int y_max):
                t_step(t_step_), x_step(x_step_), y_step(y_step_) {
            size_t t_size = ceil(double(t_max) / t_step);
            size_t x_size = ceil(double(x_max) / x_step);
            size_t y_size = ceil(double(y_max) / y_step);
            completed.resize(t_size, x_size, y_size);

            for (int current_t_start = 0; current_t_start < t_max; current_t_start += t_step) {
                int current_t_end;
                if ((current_t_start + t_step) <= t_max)
                    current_t_end = current_t_start + t_step;
                else
                    current_t_end = t_max;

                for (int current_x_start = 0; current_x_start < x_max; current_x_start += x_step) {
                    int current_x_end;
                    if ((current_x_start + x_step) <= x_max)
                        current_x_end = current_x_start + x_step;
                    else
                        current_x_end = x_max;

                    for (int current_y_start = 0; current_y_start < y_max; current_y_start += y_step) {
                        int current_y_end;
                        if ((current_y_start + y_step) <= y_max)
                            current_y_end = current_y_start + y_step;
                        else
                            current_y_end = y_max;

                        ZetaGenerationBlock current_block;
                        current_block.t_start = current_t_start;
                        current_block.t_end = current_t_end;
                        current_block.x_start = current_x_start;
                        current_block.x_end = current_x_end;
                        current_block.y_start = current_y_start;
                        current_block.y_end = current_y_end;
                        current_block.t_id = current_t_start / t_step;
                        current_block.x_id = current_x_start / x_step;
                        current_block.y_id = current_y_start / y_step;

                        queue.push_back(current_block);
                    }
                }
            }
        }

        ZetaGenerationController(ZetaGenerationController && other):
                completed(move(other.completed)), queue(move(other.queue)), t_step(other.t_step), x_step(other.x_step),
                y_step(other.y_step), mtx(move(other.mtx)){}

        bool find_available(ZetaGenerationBlock &available_block) {
            lock_guard <mutex> lock(mtx);

            if (queue.size() == 0) {
                return false;
            }

            for (auto current = queue.begin(); current != queue.end(); ++current) {
                if (is_available(*current)) {
                    available_block = *current;
                    queue.erase(current);
                    return true;
                }
            }

            return false;
        }

        bool is_available(ZetaGenerationBlock &block) {
            int t_id = block.t_id;
            int x_id = block.x_id;
            int y_id = block.y_id;

            int t_id_prev = block.t_id - 1;
            int x_id_prev = block.x_id - 1;
            int y_id_prev = block.y_id - 1;

            lock_guard <mutex> lock(mtx);

            if (t_id_prev >= 0 && !completed(t_id_prev, x_id, y_id)) {
                return false;
            }

            if (x_id_prev >= 0 && !completed(t_id, x_id_prev, y_id)) {
                return false;
            }

            if (y_id_prev >= 0 && !completed(t_id, x_id, y_id_prev)) {
                return false;
            }

            if (t_id_prev >= 0 && x_id_prev >= 0 && !completed(t_id_prev, x_id_prev, y_id)) {
                return false;
            }

            if (t_id_prev >= 0 && y_id_prev >= 0 && !completed(t_id_prev, x_id, y_id_prev)) {
                return false;
            }

            if (x_id_prev >= 0 && y_id_prev >= 0 && !completed(t_id, x_id_prev, y_id_prev)) {
                return false;
            }

            if (t_id_prev >= 0 && x_id_prev >= 0 && y_id_prev >= 0 &&
                !completed(t_id_prev, x_id_prev, y_id_prev)) {
                return false;
            }

            return true;
        }

        void set_completed(ZetaGenerationBlock &block) {
            lock_guard <mutex> lock(mtx);
            completed(block.t_id, block.x_id, block.y_id) = true;
        }

        blitz::Array<bool, 3> completed;
        list <ZetaGenerationBlock> queue;

        int t_step;
        int x_step;
        int y_step;

        mutex mtx;
    };
}
#endif //HPC_GENERATE_ZETA_PARALLEL_HH
