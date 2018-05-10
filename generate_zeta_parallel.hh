#ifndef HPC_GENERATE_ZETA_PARALLEL_HH
#define HPC_GENERATE_ZETA_PARALLEL_HH

#include <vector>
#include <mutex>
#include <unordered_set>
#include <list>

namespace parallel
{
    using namespace std;

    struct ZetaGenerationBlock
    {
        static unordered_set<tuple<int, int, int>> completed;
        static list<ZetaGenerationBlock> queue;

        static int t_step;
        static int x_step;
        static int y_step;

        static mutex mtx;

        static bool find_available(ZetaGenerationBlock& available_block){
            lock_guard<mutex> lock(mtx);

            if (queue.size() == 0){
                return false;
            }

            for(auto current = queue.begin(); current != queue.end(); ++current){
                if (current->is_available()){
                    available_block = *current;
                    queue.erase(current);
                    return true;
                }
            }

            return false;
        }

        static void initialize(int t_step, int x_step, int y_step, int t_max, int x_max, int y_max){
            ZetaGenerationBlock::t_step = t_step;
            ZetaGenerationBlock::x_step = x_step;
            ZetaGenerationBlock::y_step = y_step;

            for (int current_t_start = 0; current_t_start < t_max; current_t_start += t_step){
                int current_t_end;
                if ((current_t_start + t_step) <= t_max)
                    current_t_end = current_t_start + t_step;
                else
                    current_t_end = t_max;

                for (int current_x_start = 0; current_x_start < x_max; current_x_start += x_step){
                    int current_x_end;
                    if ((current_x_start + x_step) <= x_max)
                        current_x_end = current_x_start + x_step;
                    else
                        current_x_end = x_max;

                    for (int current_y_start = 0; current_y_start < y_max; current_y_start += y_step){
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

                        queue.push_back(current_block);
                    }
                }
            }
        }

        bool is_available(){
            int t_start_prev = t_start - t_step;
            int x_start_prev = x_start - x_step;
            int y_start_prev = y_start - y_step;

            lock_guard<mutex> lock(mtx);

            if (t_start_prev >= 0 && completed.find(make_tuple(t_start_prev, x_start, y_start)) == completed.end()){
                return false;
            }

            if (x_start_prev >= 0 && completed.find(make_tuple(t_start, x_start_prev, y_start)) == completed.end()){
                return false;
            }

            if (y_start_prev >= 0 && completed.find(make_tuple(t_start, x_start, y_start_prev)) == completed.end()){
                return false;
            }

            if (t_start_prev >= 0 && x_start_prev >= 0 &&
                completed.find(make_tuple(t_start_prev, x_start_prev, y_start)) == completed.end()){
                return false;
            }

            if (t_start_prev >= 0 && y_start_prev >= 0 &&
                completed.find(make_tuple(t_start_prev, x_start, y_start_prev)) == completed.end()){
                return false;
            }

            if (x_start_prev >= 0 && y_start_prev >= 0 &&
                completed.find(make_tuple(t_start, x_start_prev, y_start_prev)) == completed.end()){
                return false;
            }

            if (t_start_prev >= 0 && x_start_prev >= 0 && y_start_prev >= 0 &&
                completed.find(make_tuple(t_start_prev, x_start_prev, y_start_prev)) == completed.end()){
                return false;
            }

            return true;
        }

        void set_completed(){
            lock_guard<mutex> lock(mtx);
            completed.insert(make_tuple(t_start, x_start, y_start));
        }

        int t_start;
        int x_start;
        int y_start;

        int t_end;
        int x_end;
        int y_end;
    };
}
#endif //HPC_GENERATE_ZETA_PARALLEL_HH
