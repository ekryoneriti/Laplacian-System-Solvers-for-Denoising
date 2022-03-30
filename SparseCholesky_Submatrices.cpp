//============================================================================
// Name        : SparseCholesky_Submatrices.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/odeint/util/ublas_wrapper.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/generator_iterator.hpp>


#include <exception>
#include <math.h>
#include <random>
#include <chrono>
#include <iostream>
#include <fstream>
#include "SparseCholesky.h"
#define BOOST_STACKTRACE_GNU_SOURCE_NOT_REQUIRED
#include <boost/stacktrace.hpp>

using namespace boost::numeric::ublas;

template <typename T>
void power(vector<T>* vec, unsigned p){
	for(int i = 0; i < vec->size(); i++){
        if((*vec)(i) != 0){
            (*vec)(i) = pow((*vec)(i), p);
        }
	}
}

template <typename T>
void update_diagonal(compressed_matrix<T>* D, compressed_vector<T>* vals) {
//     typedef boost::numeric::ublas::compressed_vector<T>::iterator it_t;
//     for(it_t it = vals->begin(); it != D->end(); it1++){
//         (*D)(count++) = *it;
//     }
	for(int i = 0; i < vals->size(); i++){
        if((*vals)(i) != 0){
            (*D)(i,i) = (*vals)(i);
        }
	}
}

template <typename T>
unsigned int nnz(vector<T>* vec){
	unsigned int count = 0;
//     for(typename compressed_vector<T>::iterator it = vec->begin(); it != vec->end(); it++){
//         count++;
//     }
	for(int i = 0; i < vec->size(); i++){
		if((*vec)(i) != 0){
			count++;
		}
	}
	return count;
}

template <typename T>
void find_nz_indices(vector<T>* vec, vector<unsigned int>* indices){
	unsigned int count = 0;
//     for(typename compressed_vector<T>::iterator it = vec->begin(); it != vec->end(); it++){
//         std::cout << "deterministic_clique_sample::find_nz_indices: " << *it << std::endl;
//         (*indices)(count++) = *it;
//         std::cout << "deterministic_clique_sample::(*indices)(count++): " << (*indices)(count) << std::endl;
//     }
	for(int i = 0; i < vec->size(); i++){
		if((*vec)(i) != 0){
			(*indices)(count++) = i;
		}
	}
}

template <typename T>
void inline copy_values(vector<T>* popul, vector<T>* weights, int rho){
	
//     for(typename compressed_vector<T>::iterator it = vec->begin(); it != vec->end(); it++){
//         std::cout << "deterministic_clique_sample::find_nz_indices: " << *it << std::endl;
//         (*indices)(count++) = *it;
//         std::cout << "deterministic_clique_sample::(*indices)(count++): " << (*indices)(count) << std::endl;
//     }
	for(int i = 0; i < weights->size(); i++){
		popul->insert(popul->end(), rho, (*weights)(i));
        std::cout << "sparse_clique_sample::copy_weights: " << *popul << std::endl;
	}
}

template <typename T>
void inline create_probs(vector<T>* l, vector<T>* r, int denominator, int copies){
	
//     for(typename compressed_vector<T>::iterator it = vec->begin(); it != vec->end(); it++){
//         std::cout << "deterministic_clique_sample::find_nz_indices: " << *it << std::endl;
//         (*indices)(count++) = *it;
//         std::cout << "deterministic_clique_sample::(*indices)(count++): " << (*indices)(count) << std::endl;
//     }
	for(int i = 0; i < r->size(); i++){
        if((*r)(i) != 0){
            l->insert(l->end(), copies, abs((*r)(i))/denominator);
        }
	}
    std::cout << "sparse_clique_sample::l: " << *l << std::endl;
}

template <typename T>
void inline deterministic_clique_sample(compressed_matrix<T>* S, unsigned int k, compressed_matrix<T>* C) {   
//     std::ofstream out("deterministic_clique_sample_out.txt");
//     std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
//     std::cout.rdbuf(out.rdbuf());
//	std::cout << "deterministic_clique_sample::start" << std::endl;
	matrix_row< compressed_matrix<T> >* rc = new matrix_row< compressed_matrix<T> >(*S, k);
	unsigned dim = S->size1();
	vector<T> w (subrange(*rc, k, rc->size()));
	w(k) = 0;
    delete rc;
	std::cout << "deterministic_clique_sample::w: " << w.size() << std::endl;
//	std::cout << "deterministic_clique_sample::compute number of neighbors" << std::endl;
//     auto start = std::chrono::steady_clock::now();
	unsigned int neighbours_num = nnz(&w);
//     auto end = std::chrono::steady_clock::now();
//     std::cout << "Elapsed time in milliseconds deterministic clique nnz: " 
// 		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
// 		<< " milliseconds" << std::endl;
	std::cout << "deterministic_clique_sample::neighbours_num: " << neighbours_num << std::endl;
//	std::cout << "deterministic_clique_sample::compute indices of neighbors" << std::endl;
	vector<unsigned int>* indices = new vector<unsigned int>(neighbours_num);
	find_nz_indices(&w, indices);
	std::cout << "deterministic_clique_sample::indices: " << indices->size() << std::endl;
// 	matrix<T> sumY = matrix<T>(dim, dim);
//	std::cout << "deterministic_clique_sample::start for loop to create clique" << std::endl;
//	std::cout << "deterministic_clique_sample::SumY: " << sumY << std::endl;
    int index_i = -1;
    int index_j = -1;
    std::cout << "deterministic_clique_sample:: index_i:" << index_i << std::endl;
    std::cout << "deterministic_clique_sample::C size: " << C->size1() << ", " << C->size2() << std::endl;
	for(int i = 0; i < indices->size(); i++){
//		std::cout << "deterministic_clique_sample::index " << i << ": " << (*indices)(i) << std::endl;
        index_i = (*indices)(i);
		compressed_vector<T> e1 (w.size(),1);
        e1(index_i) = 1;
		T we1 = (*S)(k,k + index_i);
// 		std::cout << "deterministic_clique_sample::found first we1: " << we1 << std::endl;
//         auto start_loop = std::chrono::steady_clock::now();
        std::cout << "deterministic_clique_sample::start inner loop for i: " << i << std::endl;
		for(int j = i+1; j < indices->size(); j++){
            std::cout << "deterministic_clique_sample::start inner loop for j: " << j << std::endl;
// 			compressed_vector<T> e2 (S->size1(),1);
//             e2((*indices)(j)) = 1;
            index_j = (*indices)(j);
			T we2 = (*S)(k,k +index_j);
			std::cout << "deterministic_clique_sample::found first we2: " << we2 << std::endl;
// 			compressed_vector<T> b12 = e1 - e2;
//			std::cout << "deterministic_clique_sample::b12: " << b12 << std::endl;
// 			matrix<T> out_prod = outer_prod (b12, b12);
//			std::cout << "deterministic_clique_sample::out_prod: " << out_prod << std::endl;
// 			T scalar = (we1 * we2 / (*S)(k,k));
//			std::cout << "deterministic_clique_sample::scalar: " << scalar << std::endl;
// 			out_prod *= scalar;
// 			std::cout << "deterministic_clique_sample::outer_prod (b12, b12): " << outer_prod (b12, b12) << std::endl;
// 			(*C) = (*C) + out_prod;
            T scalar = (we1 * we2 / (*S)(k,k));
//             std::cout << "deterministic_clique_sample:scalar: " << scalar << std::endl;
//             std::cout << "deterministic_clique_sample:(*indices)(i): " << (*indices)(i) << std::endl;
//             std::cout << "deterministic_clique_sample:(*indices)(j): " << (*indices)(j) << std::endl;
//             auto start = std::chrono::high_resolution_clock::now();
            (*C)(index_i,index_i) = (*C)(index_i,index_i) + scalar;
            (*C)(index_j,index_j) = (*C)(index_j,index_j) + scalar;
            (*C)(index_i,index_j) = (*C)(index_i,index_j) - scalar;
            (*C)(index_j,index_i) = (*C)(index_j,index_i) - scalar;
//             auto end = std::chrono::high_resolution_clock::now();
//             std::cout << "Elapsed time in milliseconds deterministic clique assignment: " 
//                 << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
//                 << " microseconds" << std::endl;
//             (*C) += (outer_prod (b12, b12) * (we1 * we2 / (*S)(k,k)));
// 			std::cout << "deterministic_clique_sample:C: " << (*C) << std::endl;
		}
//         auto end_loop = std::chrono::high_resolution_clock::now();
//         std::cout << "Elapsed time in microseconds deterministic clique for loop: " 
//             << std::chrono::duration_cast<std::chrono::microseconds>(end_loop - start_loop).count()
//             << " microseconds" << std::endl;
// 		matrix<T> out_prod2 = outer_prod (e1, e1);
//		std::cout << "deterministic_clique_sample::e1: " << e1 << std::endl;
//		std::cout << "deterministic_clique_sample::out_prod2: " << out_prod2 << std::endl;
//		std::cout << "deterministic_clique_sample::w((*indices)(i)): " << w((*indices)(i)) << std::endl;
//		std::cout << "deterministic_clique_sample::abs(w(i)) - w(i): " << abs(w((*indices)(i))) - w((*indices)(i)) << std::endl;
//		std::cout << "deterministic_clique_sample::pow(abs(w(i)) - w(i), 2): " << abs(w((*indices)(i))) - pow(w((*indices)(i)), 2) << std::endl;
//		std::cout << "deterministic_clique_sample::(pow(abs(w(i)) - w(i), 2)/(*S)(k,k)): " << (abs(w((*indices)(i))) - pow(w((*indices)(i)), 2)/(*S)(k,k)) << std::endl;
//		std::cout << "deterministic_clique_sample::sumY((*indices)(i),(*indices)(i)): " << sumY((*indices)(i),(*indices)(i)) << std::endl;
// 		T diff = (abs(w((*indices)(i))) - pow(w((*indices)(i)), 2)/(*S)(k,k) - (*C)((*indices)(i),(*indices)(i)));
//		std::cout << "deterministic_clique_sample::diff: " << diff << std::endl;
//		std::cout << "deterministic_clique_sample::diff * out_prod2: " << diff * out_prod2 << std::endl;
// 		(*C) +=  diff * out_prod2;
//         (*C) +=  (abs(w((*indices)(i))) - pow(w((*indices)(i)), 2)/(*S)(k,k) - (*C)((*indices)(i),(*indices)(i))) * outer_prod (e1, e1);
        std::cout << "deterministic_clique_sample::update clique" << std::endl;
        (*C)(index_i,index_i) =  (abs(w(index_i)) - pow(w(index_i), 2)/(*S)(k,k));
//         std::cout << "deterministic_clique_sample::C: " << (*C) << std::endl;
	}
    std::cout << "deterministic_clique_sample::delete indices" << std::endl;
    delete indices; 
    std::cout << "deterministic_clique_sample::indices deleted" << std::endl;
//     //reset to standard input again
//     std::cout.rdbuf(coutbuf);
// 	(*C) = sumY;
	std::cout << "deterministic_clique_sample::finished" << std::endl;
}

template <typename T>
void deterministic_cholesky( compressed_matrix<T>* S, compressed_matrix<T>* L, compressed_matrix<T>* D) {
    std::ofstream out("deterministic_cholesky_out.txt");
    std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
    std::cout.rdbuf(out.rdbuf());
//	std::cout << "deterministic_cholesky::start" << std::endl;
	int dim = S->size1(); // the matrix is symmetric
// 	compressed_matrix<T>* S = new compressed_matrix<T>(*A);
//     std::cout << "A: " << (*A) << std::endl;
// 	std::cout << "deterministic_cholesky::for loop in deterministic_cholesky: dim - 1: " << dim - 1 << std::endl;
    compressed_matrix<T>* C;
    compressed_matrix<T>* diag;
    compressed_matrix<T>* Lk;
    matrix_column<compressed_matrix<T> >* a;
    matrix_row<compressed_matrix<T> >* rlk;
    matrix_column<compressed_matrix<T> >* clk;
    matrix_column<compressed_matrix<T> >* lc;
    matrix_column<compressed_matrix<T> >* sc;
	for(int k = 0; k < dim - 1; k++)
	{
//         std::cout << "deterministic_cholesky::k " << k << std::endl;
		(*D)(k,k) = (*S)(k,k);
//		std::cout << "deterministic_cholesky::value of diagonal " << (*D)(k,k) << std::endl;
		if((*D)(k,k) != 0)
		{
			lc = new matrix_column<compressed_matrix<T> >(*L, k);
//			std::cout << "deterministic_cholesky::before update lc: " << lc << std::endl;
			sc = new matrix_column<compressed_matrix<T> >(*S, k);
//			std::cout << "deterministic_cholesky::before update sc: " << sc << std::endl;
			(*lc) = (*sc) / (*D)(k,k);
//			std::cout << "deterministic_cholesky::after update lc: " << lc << std::endl;
//			std::cout << "deterministic_cholesky::L: " << (*L) << std::endl;
            delete lc;
            delete sc;
		}
		C = new compressed_matrix<T>(dim-k, dim-k);
		std::cout << "deterministic_cholesky::Run deterministic_clique_sample for k: " << k << std::endl;
        auto start = std::chrono::steady_clock::now();
        std::cout << "deterministic_cholesky::Size of C: " << C->size1() << ", " << C->size2() << std::endl;
		deterministic_clique_sample(S, k, C);
        std::cout << "deterministic_cholesky::Finished running deterministic_clique_sample for k: " << k << std::endl;
        auto end = std::chrono::steady_clock::now();
        std::cout << "Elapsed time in milliseconds deterministic_clique_sample: " 
		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
		<< " milseconds" << std::endl;
// 		std::cout << "deterministic_cholesky::C: " << (*C) << std::endl;
//		std::cout << "deterministic_cholesky::compute Lk" << std::endl;
//         auto start = std::chrono::steady_clock::now();
		a = new matrix_column<compressed_matrix<T> >(*S, k);
// 		std::cout << "deterministic_cholesky::column of Lk: " << a << std::endl;
		compressed_vector<T> s (subrange(*a, k, dim));
// 		std::cout << "deterministic_cholesky::update s: " << s << std::endl;
		s *= -1;
		s(0) = 0;
		diag = new compressed_matrix<T>(dim-k, dim-k);
//		std::cout << "deterministic_cholesky::update_diagonal: " << s << std::endl;
		update_diagonal(diag, &s);
//		std::cout << "deterministic_cholesky::diagonal: " << (*diag) << std::endl;
		Lk = new compressed_matrix<T>(dim-k, dim-k);
//		std::cout << "deterministic_cholesky::append diagonal to Lk: " << Lk->size1() << " Lk size2: " << Lk->size2() << std::endl;
//		std::cout << "deterministic_cholesky::append diagonal to Lk: " << diag->size1() << " diag size2: " << diag->size2() << std::endl;
		rlk = new matrix_row<compressed_matrix<T> >(*Lk, k);
		clk = new matrix_column<compressed_matrix<T> >(*Lk, k);
		(*rlk) = -s;
		(*clk) = -s;
		(*Lk) += (*diag);
		(*Lk)(k,k) = (*D)(k,k);
//         auto end = std::chrono::steady_clock::now();
//         std::cout << "Elapsed time in milliseconds deterministic clique for Lk: " 
//             << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
//             << " milliseconds" << std::endl;
//		std::cout << "deterministic_cholesky::get mc and rc " << std::endl;
//		std::cout << "deterministic_cholesky::Lk: " << (*Lk) << std::endl;
// 		std::cout << "deterministic_cholesky::Update S: " << (*S) << std::endl;
// 		std::cout << "deterministic_cholesky::Update C: " << (*C) << std::endl;
// 		std::cout << "deterministic_cholesky::Update Lk: " << (*Lk) << std::endl;
// 		std::cout << "deterministic_cholesky::- (*Lk) + (*C): " <<  (*C) + (-1)*(*Lk) << std::endl;
        std::cout << "deterministic_cholesky::size Lk: " << Lk->size1() << ", " << Lk->size2() << std::endl;
        std::cout << "deterministic_cholesky::size S: " << C->size1() << ", " << C->size2() << std::endl;
        std::cout << "deterministic_cholesky::size S: " << S->size1() << ", " << S->size2() << std::endl;
        std::cout << "deterministic_cholesky::k: " << k << ", S->size1(): " << S->size1() << std::endl;
        std::cout << "deterministic_cholesky::k: " << k << ", S->size2(): " << S->size2()<< std::endl;
		subrange(*S, k, S->size1(), k, S->size2()) += -(*Lk);
        subrange(*S, k, S->size1(), k, S->size2()) += (*C);
        Lk->resize(0,0,false);
        C->resize(0,0,false);
        delete Lk;
        delete C;
        delete a;
        delete rlk;
        delete clk;
        delete diag;
// 		std::cout << "deterministic_cholesky::S: " << (*S) << std::endl;
	}
	std::cout << "deterministic_cholesky::update final values" << std::endl;
	(*D)(dim-1,dim-1) = (*S)(dim-1,dim-1);
	(*L)(dim-1,dim-1) = 1;
    //reset to standard input again
    std::cout.rdbuf(coutbuf);
	std::cout << "deterministic_cholesky::D: " << (*D) << std::endl;
	std::cout << "deterministic_cholesky::L: " << (*L) << std::endl;
}

template <typename T>
void inline sparse_clique_sample(compressed_matrix<T>* S, unsigned int k, compressed_matrix<T>* C, int rho) {
//	std::cout << "sparse_clique_sample::start" << std::endl;
	matrix_row< compressed_matrix<T> >* rc = new matrix_row< compressed_matrix<T> >(*S, k);
	unsigned dim = S->size1();
	vector<T> w (*rc);
	w(k) = 0;
// 	std::cout << "sparse_clique_sample::w: " << w << std::endl;
//	std::cout << "sparse_clique_sample::compute number of neighbors" << std::endl;
	unsigned int neighbours_num = nnz(&w);
//	std::cout << "sparse_clique_sample::neighbours_num: " << neighbours_num << std::endl;
//	std::cout << "sparse_clique_sample::compute indices of neighbors" << std::endl;
	vector<unsigned int>* indices = new vector<unsigned int>(neighbours_num);
	find_nz_indices(&w, indices);
    // create rho copies of the weights
    vector<T>* popul = new vector<T>();
    copy_values(&indices, &popul, rho); // population 
    T weight = accumulate(w.begin(),w.end(),0);// wsv
    std::cout << "sparse_clique_sample::weight: " << weight << std::endl;
    vector<T>* p = new vector<T>();
    create_probs(p, &w, rho*weight, rho);
    // define discrete distribution
    typedef boost::minstd_rand base_generator_type;
    base_generator_type generator(12);
    typedef boost::random::discrete_distribution<> distribution_type;
    typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;
    gen_type nonuni_gen(generator, distribution_type(p));
//     boost::generator_iterator<gen_type> nonuni_die(&nonuni_gen);
//     boost::random::random_device rng;
//     std::discrete_distribution<> d(p);
//     typedef boost::variate_generator<base_generator_type&, std::discrete_distribution<>> non_unif_gen_type;
//     non_unif_gen_type nonuni_gen(generator, distribution_type(1, neighbours_num));
//     boost::generator_iterator<gen_type> nonuni_die(&nonuni_gen);
//     non_unif_gen_type prob_gen(nonuni_gen);
    // define uniform random distribution
    
    typedef boost::uniform_int<> uniform_distribution_type;
    typedef boost::variate_generator<base_generator_type&, uniform_distribution_type> uniform_gen_type;
    uniform_gen_type uni_gen(generator, uniform_distribution_type(1, popul->size()));
//     boost::generator_iterator<uniform_gen_type> uni_die(&uni_gen);
    int e1 = -1;
    int e2 = -1;
    T weight1, weight2;
    T scalar;
    for(int i = 0; i < neighbours_num; i++){
        e1 = (*popul(nonuni_gen())); // index of 
        e2 = (*popul(uni_gen()));
        std::cout << "sparse_clique_sample::e1: " << e1 << std::endl;
        std::cout << "sparse_clique_sample::e2: " << e2 << std::endl;
        if(e1 != e2){
            weight1 = S(k,e1)/rho;
            weight2 = S(k,e2)/rho;
            scalar = (weight1*weight2)/(weight1+weight2);
            std::cout << "sparse_clique_sample::C before update: " << (*C) << std::endl;
            (*C)(e1,e1) += scalar;
            (*C)(e2,e2) += scalar;
            (*C)(e1,e2) -= scalar;
            (*C)(e1,e2) -= scalar;
            std::cout << "sparse_clique_sample::C after update: " << (*C) << std::endl;
        }
    }
    std::cout << "sparse_clique_sample::C: " << (*C) << std::endl;
// 	std::cout << "sparse_clique_sample::indices: " << (*indices) << std::endl;
// 	matrix<T> sumY = matrix<T>(dim, dim);
//	std::cout << "sparse_clique_sample::start for loop to create clique" << std::endl;
//	std::cout << "sparse_clique_sample::SumY: " << sumY << std::endl;
	for(int i = 0; i < indices->size(); i++){
//		std::cout << "sparse_clique_sample::index " << i << ": " << (*indices)(i) << std::endl;
// 		matrix<T> out_prod2 = outer_prod (e1, e1);
//		std::cout << "sparse_clique_sample::e1: " << e1 << std::endl;
//		std::cout << "sparse_clique_sample::out_prod2: " << out_prod2 << std::endl;
//		std::cout << "sparse_clique_sample::w((*indices)(i)): " << w((*indices)(i)) << std::endl;
//		std::cout << "sparse_clique_sample::abs(w(i)) - w(i): " << abs(w((*indices)(i))) - w((*indices)(i)) << std::endl;
//		std::cout << "sparse_clique_sample::pow(abs(w(i)) - w(i), 2): " << abs(w((*indices)(i))) - pow(w((*indices)(i)), 2) << std::endl;
//		std::cout << "sparse_clique_sample::(pow(abs(w(i)) - w(i), 2)/(*S)(k,k)): " << (abs(w((*indices)(i))) - pow(w((*indices)(i)), 2)/(*S)(k,k)) << std::endl;
//		std::cout << "sparse_clique_sample::sumY((*indices)(i),(*indices)(i)): " << sumY((*indices)(i),(*indices)(i)) << std::endl;
// 		T diff = (abs(w((*indices)(i))) - pow(w((*indices)(i)), 2)/(*S)(k,k) - (*C)((*indices)(i),(*indices)(i)));
//		std::cout << "sparse_clique_sample::diff: " << diff << std::endl;
//		std::cout << "sparse_clique_sample::diff * out_prod2: " << diff * out_prod2 << std::endl;
// 		(*C) +=  diff * out_prod2;
//         (*C) +=  (abs(w((*indices)(i))) - pow(w((*indices)(i)), 2)/(*S)(k,k) - (*C)((*indices)(i),(*indices)(i))) * outer_prod (e1, e1);
        (*C)((*indices)(i),(*indices)(i)) +=  (abs(w((*indices)(i))) - pow(w((*indices)(i)), 2)/(*S)(k,k) - (*C)((*indices)(i),(*indices)(i)));
        std::cout << "sparse_clique_sample::C: " << (*C) << std::endl;
	}
// 	(*C) = sumY;
//	std::cout << "sparse_clique_sample::finished" << std::endl;
}

template <typename T>
void test(compressed_matrix<T>* S, compressed_matrix<T>* L, compressed_matrix<T>* D) {
//	std::cout << "deterministic_cholesky::start" << std::endl;
	int dim = S->size1(); // the matrix is symmetric
// 	compressed_matrix<T>* S = new compressed_matrix<T>(*A);
//     std::cout << "A: " << (*A) << std::endl;
// 	std::cout << "deterministic_cholesky::for loop in deterministic_cholesky: dim - 1: " << dim - 1 << std::endl;
    compressed_matrix<T>* C;
    compressed_matrix<T>* diag;
    compressed_matrix<T>* Lk;
    matrix_column<compressed_matrix<T> >* a;
    matrix_row<compressed_matrix<T> >* rlk;
    matrix_column<compressed_matrix<T> >* clk;
    matrix_column<compressed_matrix<T> >* lc;
    matrix_column<compressed_matrix<T> >* sc;
    int rho = 1000;
	for(int k = 0; k < dim - 1; k++)
	{
//         std::cout << "deterministic_cholesky::k " << k << std::endl;
		(*D)(k,k) = (*S)(k,k);
//		std::cout << "deterministic_cholesky::value of diagonal " << (*D)(k,k) << std::endl;
		if((*D)(k,k) != 0)
		{
			lc = new matrix_column<compressed_matrix<T> >(*L, k);
//			std::cout << "deterministic_cholesky::before update lc: " << lc << std::endl;
			sc = new matrix_column<compressed_matrix<T> >(*S, k);
//			std::cout << "deterministic_cholesky::before update sc: " << sc << std::endl;
			(*lc) = (*sc) / (*D)(k,k);
//			std::cout << "deterministic_cholesky::after update lc: " << lc << std::endl;
//			std::cout << "deterministic_cholesky::L: " << (*L) << std::endl;
		}
		C = new compressed_matrix<T>(dim, dim,k);
//		std::cout << "deterministic_cholesky::Run deterministic_clique_sample" << std::endl;
//         auto start = std::chrono::steady_clock::now();
       
		sparse_clique_sample(S, k, C, rho);
//         auto end = std::chrono::steady_clock::now();
        
//         std::cout << "Elapsed time in milliseconds : " 
// 		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
// 		<< " milseconds" << std::endl;
// 		std::cout << "deterministic_cholesky::C: " << (*C) << std::endl;
//		std::cout << "deterministic_cholesky::compute Lk" << std::endl;
		a = new matrix_column<compressed_matrix<T> >(*S, k);
//		std::cout << "deterministic_cholesky::column of Lk: " << a << std::endl;
		compressed_vector<T> s (*a);
//		std::cout << "deterministic_cholesky::update s: " << a << std::endl;
		s *= -1;
		s(k) = 0;
		diag = new compressed_matrix<T>(dim, dim);
//		std::cout << "deterministic_cholesky::update_diagonal: " << s << std::endl;
		update_diagonal(diag, &s);
//		std::cout << "deterministic_cholesky::diagonal: " << (*diag) << std::endl;
		Lk = new compressed_matrix<T>(dim, dim);
//		std::cout << "deterministic_cholesky::append diagonal to Lk: " << Lk->size1() << " Lk size2: " << Lk->size2() << std::endl;
//		std::cout << "deterministic_cholesky::append diagonal to Lk: " << diag->size1() << " diag size2: " << diag->size2() << std::endl;
		rlk = new matrix_row<compressed_matrix<T> >(*Lk, k);
		clk = new matrix_column<compressed_matrix<T> >(*Lk, k);
		(*rlk) = -s;
		(*clk) = -s;
		(*Lk) += (*diag);
		(*Lk)(k,k) = (*D)(k,k);
//		std::cout << "deterministic_cholesky::get mc and rc " << std::endl;
//		std::cout << "deterministic_cholesky::Lk: " << (*Lk) << std::endl;
//		std::cout << "deterministic_cholesky::Update S: " << (*S) << std::endl;
//		std::cout << "deterministic_cholesky::Update C: " << (*C) << std::endl;
//		std::cout << "deterministic_cholesky::Update Lk: " << (*Lk) << std::endl;
//		std::cout << "deterministic_cholesky::- (*Lk) + (*C): " <<  (*C) + (-1)*(*Lk) << std::endl;
		(*S) += -(*Lk);
        (*S) += (*C);
// 		std::cout << "deterministic_cholesky::S: " << (*S) << std::endl;
	}
//	std::cout << "deterministic_cholesky::update final values" << std::endl;
	(*D)(dim-1,dim-1) = (*S)(dim-1,dim-1);
	(*L)(dim-1,dim-1) = 1;
// 	std::cout << "deterministic_cholesky::D: " << (*D) << std::endl;
// 	std::cout << "deterministic_cholesky::L: " << (*L) << std::endl;
}

int main () {
	try{
		compressed_matrix<double> m (3,3);
		m(0,0) = 100;
		m(0,1) = -50;
		m(1,0) = -50;
		m(0,2) = -90;
		m(2,0) = -90;
		m(1,1) = 50;
		m(1,2) = 95;
		m(2,1) = 95;
		m(2,2) = 350;
		//
//		m(0,0) = 9;
//		m(0,1) = 3;
//		m(1,0) = 3;
//		m(1,1) = 5;
//		m(1,2) = 2;
//		m(2,1) = 2;
//		m(2,2) = 17;
		std::cout << "m: " << m << std::endl;
		compressed_matrix<double> L(3,3);
		compressed_matrix<double> D(3, 3);
		deterministic_cholesky(&m, &L, &D);
		std::cout << "L: " << L << std::endl;
		std::cout << "D: " << D << std::endl;
	}
	catch(std::exception& e){
		std::cout << e.what() << std::endl;
		std::cout << boost::stacktrace::stacktrace();
	}
}
