//============================================================================
// Name        : SparseCholesky_Mat.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

// #include <chrono>
#include <matrix.h>
#include <mex.h>
#include "SparseCholesky_Mat.h"

/* Definitions to keep compatibility with earlier versions of ML */
#ifndef MWSIZE_MAX
typedef int mwSize;
typedef int mwIndex;
typedef int mwSignedIndex;
#if (defined(_LP64) || defined(_WIN64)) && !defined(MX_COMPAT_32)
/* Currently 2^48 based on hardware limitations */
# define MWSIZE_MAX    281474976710655UL
# define MWINDEX_MAX   281474976710655UL
# define MWSINDEX_MAX  281474976710655L
# define MWSINDEX_MIN -281474976710655L
#else
# define MWSIZE_MAX    2147483647UL
# define MWINDEX_MAX   2147483647UL
# define MWSINDEX_MAX  2147483647L
# define MWSINDEX_MIN -2147483647L
#endif
#define MWSIZE_MIN    0UL
#define MWINDEX_MIN   0UL
#endif

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
// #include "SparseCholesky_Mat.h"
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
void sparse_cholesky( boost::numeric::ublas::compressed_matrix<T>* S, boost::numeric::ublas::compressed_matrix<T>* L, boost::numeric::ublas::compressed_matrix<T>* D) {
	std::cout << "sparse_cholesky::start" << std::endl;
	int dim = S->size1(); // the matrix is symmetric
// 	compressed_matrix<T>* S = new compressed_matrix<T>(*A);
    std::cout << "S: " << (*S) << std::endl;
// 	std::cout << "sparse_cholesky::for loop in sparse_cholesky: dim - 1: " << dim - 1 << std::endl;
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
//         std::cout << "sparse_cholesky::k " << k << std::endl;
		(*D)(k,k) = (*S)(k,k);
//		std::cout << "sparse_cholesky::value of diagonal " << (*D)(k,k) << std::endl;
		if((*D)(k,k) != 0)
		{
			lc = new matrix_column<compressed_matrix<T> >(*L, k);
//			std::cout << "sparse_cholesky::before update lc: " << lc << std::endl;
			sc = new matrix_column<compressed_matrix<T> >(*S, k);
//			std::cout << "sparse_cholesky::before update sc: " << sc << std::endl;
			(*lc) = (*sc) / (*D)(k,k);
//			std::cout << "sparse_cholesky::after update lc: " << lc << std::endl;
//			std::cout << "sparse_cholesky::L: " << (*L) << std::endl;
		}
		C = new compressed_matrix<T>(dim, dim,k);
		std::cout << "sparse_cholesky::Run sparse_clique_sample" << std::endl;
//         auto start = std::chrono::steady_clock::now();
       
		sparse_clique_sample(S, k, C, rho);
//         auto end = std::chrono::steady_clock::now();
        
//         std::cout << "Elapsed time in milliseconds : " 
// 		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
// 		<< " milseconds" << std::endl;
// 		std::cout << "sparse_cholesky::C: " << (*C) << std::endl;
//		std::cout << "sparse_cholesky::compute Lk" << std::endl;
		a = new matrix_column<compressed_matrix<T> >(*S, k);
//		std::cout << "sparse_cholesky::column of Lk: " << a << std::endl;
		compressed_vector<T> s (*a);
//		std::cout << "sparse_cholesky::update s: " << a << std::endl;
		s *= -1;
		s(k) = 0;
		diag = new compressed_matrix<T>(dim, dim);
//		std::cout << "sparse_cholesky::update_diagonal: " << s << std::endl;
		update_diagonal(diag, &s);
//		std::cout << "sparse_cholesky::diagonal: " << (*diag) << std::endl;
		Lk = new compressed_matrix<T>(dim, dim);
//		std::cout << "sparse_cholesky::append diagonal to Lk: " << Lk->size1() << " Lk size2: " << Lk->size2() << std::endl;
//		std::cout << "sparse_cholesky::append diagonal to Lk: " << diag->size1() << " diag size2: " << diag->size2() << std::endl;
		rlk = new matrix_row<compressed_matrix<T> >(*Lk, k);
		clk = new matrix_column<compressed_matrix<T> >(*Lk, k);
		(*rlk) = -s;
		(*clk) = -s;
		(*Lk) += (*diag);
		(*Lk)(k,k) = (*D)(k,k);
//		std::cout << "sparse_cholesky::get mc and rc " << std::endl;
//		std::cout << "sparse_cholesky::Lk: " << (*Lk) << std::endl;
//		std::cout << "sparse_cholesky::Update S: " << (*S) << std::endl;
//		std::cout << "sparse_cholesky::Update C: " << (*C) << std::endl;
//		std::cout << "sparse_cholesky::Update Lk: " << (*Lk) << std::endl;
//		std::cout << "sparse_cholesky::- (*Lk) + (*C): " <<  (*C) + (-1)*(*Lk) << std::endl;
		(*S) += -(*Lk);
        (*S) += (*C);
		std::cout << "sparse_cholesky::S: " << (*S) << std::endl;
	}
//	std::cout << "sparse_cholesky::update final values" << std::endl;
	(*D)(dim-1,dim-1) = (*S)(dim-1,dim-1);
	(*L)(dim-1,dim-1) = 1;
// 	std::cout << "sparse_cholesky::D: " << (*D) << std::endl;
// 	std::cout << "sparse_cholesky::L: " << (*L) << std::endl;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    auto start3 = std::chrono::steady_clock::now();
/* check for proper number of arguments */
    if(nrhs!=5) {
        mexErrMsgIdAndTxt("sparse_cholesky_matlab_api:nrhs","Five inputs required.");
    }
    if(nlhs!=6) {
        mexErrMsgIdAndTxt("sparse_cholesky_matlab_api:nlhs","Six outputs required.");
    }
//declare variables
//     mxArray *l_out_m, *d_out_m;
    double * x_in_m;
    double * y_in_m;
    int rows, cols;
    double * a_in_m;
    const mwSize *dims;
    double *l, *d;
    int dimx, numdims;
 
//associate pointers
    x_in_m = mxGetPr(prhs[0]);
    y_in_m = mxGetPr(prhs[1]);
    a_in_m = mxGetPr(prhs[2]); // input A
    /* Get the size and pointers to input data */
    dims = mxGetDimensions(prhs[0]);
    numdims = mxGetNumberOfDimensions(prhs[0]);
    dimx = (int)dims[0];
//     std::cout << "dimx: " << dimx << std::endl;
    rows = (int) mxGetScalar(prhs[3]); 
    cols = (int) mxGetScalar(prhs[4]);
    
    using namespace boost::numeric::ublas;
    auto start2 = std::chrono::steady_clock::now();
    compressed_matrix<double>* A = new compressed_matrix<double>(rows, cols, rows);
    for(int i = 0; i < dimx; i++){
//         std::cout << i << ": " << x_in_m[i] << ", " << y_in_m[i] << ": " << a_in_m[i] << std::endl;
        (*A)(x_in_m[i]-1,y_in_m[i]-1) = a_in_m[i];
//         std::cout << (*A)(x_in_m[i]-1,y_in_m[i]-1) << std::endl;
    }
    
    auto end2 = std::chrono::steady_clock::now();
    std::cout << "Elapsed time in seconds copying data : " 
		<< std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count()
		<< " mil" << std::endl;
//     
//     l = mxGetPr(l_out_m);
//     d = mxGetPr(d_out_m);

//associate outputs
//     l_out_m = plhs[0] = mxCreateDoubleMatrix(rows,cols,mxREAL);
//     d_out_m = plhs[1] = mxCreateDoubleMatrix(rows,cols,mxREAL);
//     plhs[0] = mxCreateSparse(rows,cols,mxREAL);
//     plhs[1] = mxCreateSparse(rows,cols,mxREAL);
    
// run sparse cholesky
    compressed_matrix<double>* L = new compressed_matrix<double>(rows, cols, rows);
    compressed_matrix<double>* D = new compressed_matrix<double>(rows, cols, rows);
    auto start = std::chrono::steady_clock::now();
    sparse_cholesky(A, L, D);
//     test();
    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed time in milliseconds sparse: " 
		<< std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
		<< " milliseconds" << std::endl;
//     std::cout << "L: " << (*L) << std::endl;
//     std::cout << "D: " << (*D) << std::endl;
    // copy data to matlab matrices
    auto start1 = std::chrono::steady_clock::now();
    typedef boost::numeric::ublas::compressed_matrix<double>::iterator1 it1_t;
    typedef boost::numeric::ublas::compressed_matrix<double>::iterator2 it2_t;
//     double *dynamicData = (double*)mxMalloc(rows * cols * sizeof( double ));
//     double *pr;
//     int *ir, *jc;
    double* data = (double*)mxMalloc(rows * cols * sizeof( double ));
    int* r_indices = (int*)mxMalloc(rows * cols * sizeof( int ));
    int* c_indices = (int*)mxMalloc(rows * cols * sizeof( int ));
    int count = 0;
    int r_count = 0;
    int c_count = 0;
    for (it1_t it1 = L->begin1(); it1 != L->end1(); it1++)
    {
      
//       std::cout << "it1.index1(): " << it1.index1() << std::endl;
      for (it2_t it2 = it1.begin(); it2 != it1.end(); it2++)
      {
        count++;
//         std::cout << "it2.index1(): " << it2.index1() << std::endl;
        r_indices[r_count++] = (int)it2.index1() + 1;
        c_indices[c_count] = (int)it2.index2() + 1;
        data[c_count++] = (double)(*it2);
//         std::cout << "it2.index2(): " << it2.index2() << std::endl;
//         dynamicData[it2.index1() * rows + it2.index2()] = *it2;
      }
    }
    plhs[0] = mxCreateDoubleMatrix(count,1,mxREAL);
    mxSetData(plhs[0], data);
    plhs[1] = mxCreateNumericMatrix(count,1,mxINT32_CLASS,mxREAL);
    mxSetData(plhs[1], r_indices);
    plhs[2] = mxCreateNumericMatrix(count,1,mxINT32_CLASS,mxREAL);
    mxSetData(plhs[2], c_indices);
//     std::cout << "c_indices: " << c_indices << std::endl;
//     std::cout << "r_indices: " << r_indices << std::endl;
//     plhs[0] = mxCreateSparse(rows,cols,count,mxREAL);
//     pr = (double *)mxGetPr (plhs[0]);
//     ir = (int *)mxGetIr (plhs[0]);
//     jc = (int *)mxGetJc (plhs[0]);
//     memcpy(pr, data, count*sizeof(double));
//     memcpy(ir, r_indices, count*sizeof(int));
//     memcpy(jc, c_indices, count*sizeof(int));
//     mxSetData(plhs[0], dynamicData);
    count = 0;
    double* data2 = (double*)mxMalloc(rows * cols * sizeof( double ));
    int* r_indices2 = (int*)mxMalloc(rows * cols * sizeof( int ));
    int* c_indices2 = (int*)mxMalloc(rows * cols * sizeof( int ));
    count = 0;
    r_count = 0;
    c_count = 0;
//     plhs[1] = mxCreateSparse(rows,cols,D->size1(),mxREAL);
//     double *dynamicData2 = (double*)mxMalloc(rows * cols * sizeof( double ));
    for (it1_t it1 = D->begin1(); it1 != D->end1(); it1++)
    {
      for (it2_t it2 = it1.begin(); it2 != it1.end(); it2++)
      {
//         std::cout << "value: " << *it2 << std::endl;
        count++;
        r_indices2[r_count++] = (int)it2.index1() + 1;
        c_indices2[c_count] = (int)it2.index2() + 1;
        data2[c_count++] = (double)(*it2);
//         dynamicData2[it2.index1() * rows + it2.index2()] = *it2;
      }
    }
    plhs[3] = mxCreateDoubleMatrix(count,1,mxREAL);
    mxSetData(plhs[3], data2);
    plhs[4] = mxCreateNumericMatrix(count,1,mxINT32_CLASS,mxREAL);
    mxSetData(plhs[4], r_indices2);
    plhs[5] = mxCreateNumericMatrix(count,1,mxINT32_CLASS,mxREAL);
    mxSetData(plhs[5], c_indices2);
//     std::cout << "dynamicData2: " << *dynamicData2 << std::endl;
    
//     mxSetData(plhs[1], dynamicData2);
    
//     std::cout << "d_out_m: " << *d_out_m << std::endl;
//     std::cout << "l_out_m: " << *l_out_m << std::endl;
    auto end1 = std::chrono::steady_clock::now();
    std::cout << "Elapsed time in milliseconds : " 
		<< std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count()
		<< " mil" << std::endl;
    auto end3 = std::chrono::steady_clock::now();
    std::cout << "Elapsed time in milliseconds : " 
		<< std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count()
		<< " mil" << std::endl;
    return;
}
