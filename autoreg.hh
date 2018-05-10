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
#include <fstream>
#include <vector>
#include <thread>
#include <mutex>

#include <blitz/array.h>         // for Array, Range, shape, any

#include "sysv.hh"               // for sysv
#include "types.hh"              // for size3, ACF, AR_coefs, Zeta, Array2D
#include "voodoo.hh"             // for generate_AC_matrix
#include "parallel_mt.hh"
#include "generate_zeta_parallel.hh"

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
		//Чтение из файла и создание генераторов
		size_t n_threads = 8;
		std::string file_name = "mt_configs_file";
		std::ifstream file(file_name);

		std::vector<parallel_mt> generators;
		for (size_t i=0; i<n_threads; ++i){
			mt_config config;
			file>>config;
			generators.push_back(parallel_mt(config));
		}

		//Параллельное заполнение массива
		std::normal_distribution<T> normal(T(0), std::sqrt(variance));
		Zeta<T> eps(size);

		auto current_start = std::begin(eps);
		size_t step = size(0)*size(1)*size(2) / n_threads;
		std::vector<std::thread> threads;

		for (size_t i=0; i<n_threads-1; ++i){
			auto current_end = std::next(current_start, step);
			auto generator = std::bind(normal, generators[i]);
			std::thread current_thread(
					std::generate<decltype(current_start), decltype
					(generator)>, current_start, current_end, generator);
			threads.push_back(std::move(current_thread));
			current_start = current_end;
		}

		auto current_end = std::end(eps);
		auto generator = std::bind(normal, generators[generators.size()-1]);
		std::thread current_thread(
				std::generate<decltype(current_start), decltype
				(generator)>, current_start, current_end, generator);
		threads.push_back(std::move(current_thread));

		//Ожидание выполнения
		for(auto& current_thread : threads){
			current_thread.join();
		}

		//Проверка
		if (std::any_of(std::begin(eps), std::end(eps), &::autoreg::isnan<T>)) {
			throw std::runtime_error("white noise generator produced some NaNs");
		}
		return eps;
	}

	void generate_zeta_parallel_worker(parallel::ZetaGenerationController& controller){
		std::cout<<"hello";
//		const size3 fsize = phi.shape();
//		const size3 zsize = zeta.shape();
//
//		while (true){
//			parallel::ZetaGenerationBlock block;
//			bool was_found = controller.find_available(block);
//
//			if (!was_found){
//				std::lock_guard<std::mutex> lock(controller.mtx);
//				if (controller.queue.size() == 0){
//					break;
//				}
//			} else {
//				for (int t=block.t_start; t<block.t_end; t++) {
//					for (int x=block.x_start; x<block.x_end; x++) {
//						for (int y=block.y_start; y<block.y_end; y++) {
//							const int m1 = std::min(t+1, fsize[0]);
//							const int m2 = std::min(x+1, fsize[1]);
//							const int m3 = std::min(y+1, fsize[2]);
//							T sum = 0;
//							for (int k=0; k<m1; k++)
//								for (int i=0; i<m2; i++)
//									for (int j=0; j<m3; j++)
//										sum += phi(k, i, j)*zeta(t-k, x-i, y-j);
//							zeta(t, x, y) += sum;
//						}
//					}
//				}
//
//				controller.set_completed(block);
//			}
//		}
	}

	/// Генерация отдельных частей реализации волновой поверхности.
	template<class T>
	void generate_zeta(const AR_coefs<T>& phi, Zeta<T>& zeta) {
		const size3 fsize = phi.shape();
		const size3 zsize = zeta.shape();
		const int t1 = zsize[0];
		const int x1 = zsize[1];
		const int y1 = zsize[2];

		const int t_step = fsize[0];
		const int x_step = fsize[1];
		const int y_step = fsize[2];

		parallel::ZetaGenerationController controller(t_step, x_step, y_step, t1, x1, y1);

		int n_threads = 8;
		std::vector<std::thread> threads;

		for(int thread_id=0; thread_id< n_threads; ++thread_id){
			std::thread current_thread(generate_zeta_parallel_worker, std::ref(controller));
			threads.push_back(std::move(current_thread));
		}

		for(std::thread& current_thread : threads){
			current_thread.join();
		}

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

#endif // AUTOREG_HH
