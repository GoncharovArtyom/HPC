#ifndef AUTOREG_HH
#define AUTOREG_HH

#include <algorithm>             // for min, any_of, copy_n, for_each, generate
#include <cassert>               // for assert
#include <chrono>                // for duration, steady_clock, steady_clock...
#include <cmath>                 // for isnan
#include <cstdlib>               // for abs
#include <functional>            // for bind
#include <iostream>              // for operator<<, cerr, endl
#include <fstream>               // for ofstream
#include <random>                // for mt19937, normal_distribution
#include <stdexcept>             // for runtime_error
#include <vector>
#include <mutex>
#include <thread>
#include <queue>
#include <ctime>
#include <blitz/array.h>         // for Array, Range, shape, any
#include "parallel_mt.hh"
#include "sysv.hh"               // for sysv
#include "types.hh"              // for size3, ACF, AR_coefs, Zeta, Array2D
#include "voodoo.hh"             // for generate_AC_matrix

/// @file
/// File with subroutines for AR model, Yule-Walker equations
/// and some others.

namespace autoreg {

    template<class T>
    ACF<T>
    approx_acf(T alpha, T beta, T gamm, const Vec3<T>& delta, const size3& acf_size) {
        ACF<T> acf(acf_size);
        blitz::firstIndex t;
        blitz::secondIndex x;
        blitz::thirdIndex y;
        acf = gamm
              * blitz::exp(-alpha * (t*delta[0] + x*delta[1] + y*delta[2]))
              //	 		* blitz::cos(beta * (t*delta[0] + x*delta[1] + y*delta[2]));
              * blitz::cos(beta * t * delta[0])
              * blitz::cos(beta * x * delta[1])
              * blitz::cos(beta * y * delta[2]);
        return acf;
    }

    template<class T>
    T white_noise_variance(const AR_coefs<T>& ar_coefs, const ACF<T>& acf) {
        return acf(0,0,0) - blitz::sum(ar_coefs * acf);
    }

    template<class T>
    T ACF_variance(const ACF<T>& acf) {
        return acf(0,0,0);
    }

    /// Удаление участков разгона из реализации.
    template<class T>
    Zeta<T>
    trim_zeta(const Zeta<T>& zeta2, const size3& zsize) {
        using blitz::Range;
        using blitz::toEnd;
        size3 zsize2 = zeta2.shape();
        return zeta2(
                Range(zsize2(0) - zsize(0), toEnd),
                Range(zsize2(1) - zsize(1), toEnd),
                Range(zsize2(2) - zsize(2), toEnd)
        );
    }

    template<class T>
    bool is_stationary(AR_coefs<T>& phi) {
        return !blitz::any(blitz::abs(phi) > T(1));
    }

    template<class T>
    AR_coefs<T>
    compute_AR_coefs(const ACF<T>& acf) {
        using blitz::Range;
        using blitz::toEnd;
        const int m = acf.numElements()-1;
        Array2D<T> acm = generate_AC_matrix(acf);
        //{ std::ofstream out("acm"); out << acm; }

        /**
        eliminate the first equation and move the first column of the remaining
        matrix to the right-hand side of the system
        */
        Array1D<T> rhs(m);
        rhs = acm(Range(1, toEnd), 0);
        //{ std::ofstream out("rhs"); out << rhs; }

        // lhs is the autocovariance matrix without first
        // column and row
        Array2D<T> lhs(blitz::shape(m,m));
        lhs = acm(Range(1, toEnd), Range(1, toEnd));
        //{ std::ofstream out("lhs"); out << lhs; }

        assert(lhs.extent(0) == m);
        assert(lhs.extent(1) == m);
        assert(rhs.extent(0) == m);
        sysv<T>('U', m, 1, lhs.data(), m, rhs.data(), m);
        AR_coefs<T> phi(acf.shape());
        assert(phi.numElements() == rhs.numElements() + 1);
        phi(0,0,0) = 0;
        std::copy_n(rhs.data(), rhs.numElements(), phi.data()+1);
        //{ std::ofstream out("ar_coefs"); out << phi; }
        if (!is_stationary(phi)) {
            std::cerr << "phi.shape() = " << phi.shape() << std::endl;
            std::for_each(
                    phi.begin(),
                    phi.end(),
                    [] (T val) {
                        if (std::abs(val) > T(1)) {
                            std::cerr << val << std::endl;
                        }
                    }
            );
            throw std::runtime_error("AR process is not stationary, i.e. |phi| > 1");
        }
        return phi;
    }

    template<class T>
    bool
    isnan(T rhs) noexcept {
        return std::isnan(rhs);
    }

    /// Генерация белого шума по алгоритму Вихря Мерсенна и
    /// преобразование его к нормальному распределению по алгоритму Бокса-Мюллера.
    template<class T>
    Zeta<T>
    generate_white_noise(const size3& size, const T variance) {
        if (variance < T(0)) {
            throw std::runtime_error("variance is less than zero");
        }
//*
        int OMP_NUM_THREADS=8;
        autoreg::mt_config *arr_mt_config = new autoreg::mt_config[OMP_NUM_THREADS];
        std::fstream configFile("mt_config", std::ios_base::in);

        //for(int i=0; i<OMP_NUM_THREADS; ++i){
        configFile >> arr_mt_config[0];
        //}
        configFile.close();

        std::vector<std::function<void(int)>> vfunc;
        Zeta<T> eps(size);
        long shift = size[0]*size[1]*size[2]/OMP_NUM_THREADS;
        for(int it=0; it<OMP_NUM_THREADS; ++it){
            vfunc.push_back([&](int num){
                autoreg::parallel_mt mt_generator(arr_mt_config[0]);
#if !defined(DISABLE_RANDOM_SEED)
                mt_generator.seed(std::chrono::steady_clock::now().time_since_epoch().count();
#endif
                std::normal_distribution<T> normal(T(0), std::sqrt(variance));

                long shift2 = shift*num;
                auto start = std::next(std::begin(eps), shift2);

                if(it != OMP_NUM_THREADS-1){
                    std::generate(start, std::next(start, shift), std::bind(normal, mt_generator));
                    if(std::any_of(start, std::next(start, shift), &autoreg::isnan<T>)){
                        throw std::runtime_error("white noise generator produced some NaNs");
                    }
                }
                else{
                    std::generate(start, std::end(eps), std::bind(normal, mt_generator));
                    if(std::any_of(start, std::end(eps), &autoreg::isnan<T>)){
                        throw std::runtime_error("white noise generator produced some NaNs");
                    }
                }
                //std::clog<<"from thread"<<std::endl;
            });
        }


        std::vector<std::thread> vThread;
        for(int it = 0; it < OMP_NUM_THREADS; ++it){
            std::thread newThread(vfunc[it], it);
            vThread.push_back(move(newThread));
        }

        for(int it=0; it<OMP_NUM_THREADS; ++it){
            //std::clog<< "wait "<< it <<std::endl;

            vThread[it].join();
        }



// */

/*

		// инициализация генератора
		std::mt19937 generator;
		#if !defined(DISABLE_RANDOM_SEED)
		generator.seed(std::chrono::steady_clock::now().time_since_epoch().count());
		#endif
		std::normal_distribution<T> normal(T(0), std::sqrt(variance));

		// генерация и проверка
		Zeta<T> eps(size);
		std::generate(std::begin(eps), std::end(eps), std::bind(normal, generator));
		if (std::any_of(std::begin(eps), std::end(eps), &::autoreg::isnan<T>)) {
			throw std::runtime_error("white noise generator produced some NaNs");
		}

// */


        return eps;
    }


    struct block{
        int x,y,t;
        int x_start, x_end, y_start, y_end, t_start, t_end;
        bool operator==(const block& rhs)
        {
            return x==rhs.x && y==rhs.y && t==rhs.t;
        }
    };
    struct availability_block{
        int x,y,t;
        bool available=false;
    };
    struct completed_pushed{
        bool completed=false;
        bool pushed=false;
    };

    bool is_ready_to_queue(int t, int x, int y, blitz::Array<completed_pushed,3> &map){
        //size3 sizes = map.shape();
        //if(map(t,x,y).pushed)
        //	return false;

        if(t==0 || map(t-1,x,y).completed)
            if(x==0 || map(t,x-1,y).completed)
                if(y==0 || map(t,x,y-1).completed)
                    return true;

        return false;
    }
    void push_block(int t, int x, int y, blitz::Array<completed_pushed,3> &map){
        //qu.push(blocks(t,x,y));
        map(t,x,y).pushed=true;
//std::cout<< t<<' '<<x<<' '<<y<<" pushed";
    }

    /// Генерация отдельных частей реализации волновой поверхности.
    template<class T>
    void generate_zeta(const AR_coefs<T>& phi, Zeta<T>& zeta) {
//*
        int num_proc = 8;
        std::vector<std::function<void()>> vfunc;
        std::vector<availability_block> map_blocks;
        std::mutex q_mutex;
        std::mutex map_mutex;

        std::queue<block> q_blocks;
        const size3 fsize = phi.shape();
        const size3 zsize = zeta.shape();
        const int t1 = zsize[0];
        const int x1 = zsize[1];
        const int y1 = zsize[2];
        int xys = 4;
        int ts = 20;
        int size_block_x = x1/xys;
        int size_block_y = y1/xys;
        int size_block_t = t1/ts;

        block bl;
        availability_block a_bl;

        size3 sizess(ts, xys, xys);
        //blitz::Array<block,3> blocks(sizess);
        std::vector<block> blocks_vec;
        blitz::Array<completed_pushed,3> map_bl(sizess);

        for(int k=0; k<ts; ++k){
            for(int i=0; i<xys; ++i){
                for(int j=0; j<xys; ++j){
                    map_bl(k,i,j)=completed_pushed();
                    bl = block();
                    a_bl = availability_block();
                    bl.x=a_bl.x=i; bl.y=a_bl.y=j; bl.t=a_bl.t=k;
                    bl.x_start = i*size_block_x; bl.x_end=(i+1)*size_block_x;
                    bl.y_start = j*size_block_y; bl.y_end=(j+1)*size_block_y;
                    bl.t_start = k*size_block_t; bl.t_end=(k+1)*size_block_t;
                    if(i==xys-1){
                        bl.x_end=x1;
                    }
                    if(j==xys-1){
                        bl.y_end=y1;
                    }
                    if(k==ts-1){
                        bl.t_end=t1;
                    }
                    //blocks(k,i,j)=bl;
                    blocks_vec.push_back(bl);
                }
            }
        }
//block b = blocks_vec[0];
        std::cout<<"count blocks is "<<blocks_vec.size()<<std::endl;

        //q_blocks.push();
        push_block(0, 0, 0, map_bl);
        bool done = false;
        std::clog<<("init threads")<<std::endl;
        for(int it=0; it<num_proc; ++it){
            vfunc.push_back([&](){
                block b;
                bool init = false;
                bool pushed = false;
                int bx=-1;
                int by=-1;
                int bt=-1;
                while(blocks_vec.size()>0){
                    q_mutex.lock();{
                        if(bx!=-1 ){
                            map_bl(bt, bx, by).completed = true;
                        }
                        for(block bl : blocks_vec){
                            if(is_ready_to_queue(bl.t, bl.x, bl.y, map_bl)){
                                b = block(bl);
                                init = true;
                                auto iter = std::find(blocks_vec.begin(), blocks_vec.end(), bl);
                                blocks_vec.erase(iter);
                                break;
                            }
                        }
                        q_mutex.unlock();
                        if(init){
                            init = false;
                            for (int t=b.t_start; t<b.t_end; t++) {
                                for (int x=b.x_start; x<b.x_end; x++) {
                                    for (int y=b.y_start; y<b.y_end; y++) {
                                        const int m1 = std::min(b.t_start+1, fsize[0]);
                                        const int m2 = std::min(b.x_start+1, fsize[1]);
                                        const int m3 = std::min(b.y_start+1, fsize[2]);
                                        T sum = 0;
                                        for (int k=0; k<m1; k++)
                                            for (int i=0; i<m2; i++)
                                                for (int j=0; j<m3; j++)
                                                    sum += phi(k, i, j)*zeta(t-k, x-i, y-j);
                                        zeta(t, x, y) += sum;
                                    }
                                }
                            }
                            bx = b.x;
                            by = b.y;
                            bt = b.t;
//*
                        }
                    }
                }

            });

        }
        std::clog<<"start thread"<<std::endl;
        std::vector<std::thread> vThread;
        for(int it = 0; it < num_proc; ++it){
            std::thread newThread(vfunc[it]);
            vThread.push_back(move(newThread));
        }
        for(int it=0; it<num_proc; ++it){
//std::clog<<"wait " <<it<<std::endl;
            vThread[it].join();
        }
//*/
/*
		const size3 fsize = phi.shape();
		const size3 zsize = zeta.shape();
		const int t1 = zsize[0];
		const int x1 = zsize[1];
		const int y1 = zsize[2];
		for (int t=0; t<t1; t++) {
			for (int x=0; x<x1; x++) {
				for (int y=0; y<y1; y++) {
					const int m1 = std::min(t+1, fsize[0]);
					const int m2 = std::min(x+1, fsize[1]);
					const int m3 = std::min(y+1, fsize[2]);
					T sum = 0;
					for (int k=0; k<m1; k++)
						for (int i=0; i<m2; i++)
							for (int j=0; j<m3; j++)
								sum += phi(k, i, j)*zeta(t-k, x-i, y-j);
					zeta(t, x, y) += sum;
				}
			}
		}
//*/
    }

    template<class T, int N>
    T mean(const blitz::Array<T,N>& rhs) {
        return blitz::sum(rhs) / rhs.numElements();
    }

    template<class T, int N>
    T variance(const blitz::Array<T,N>& rhs) {
        assert(rhs.numElements() > 0);
        const T m = mean(rhs);
        return blitz::sum(blitz::pow(rhs-m, 2)) / (rhs.numElements() - 1);
    }


}
