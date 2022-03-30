/*
 * SparceCholesky.h
 *
 *  Created on: Jul. 14, 2019
 *      Author: evangelia
 */



#ifndef SPARSECHOLESKY_H_
#define SPARSECHOLESKY_H_

#include <boost/numeric/ublas/matrix_sparse.hpp>

using namespace boost::numeric::ublas;


template <typename T>
void deterministic_cholesky(compressed_matrix<T>* A, compressed_matrix<T>* L, compressed_matrix<T>* D);

template <typename T>
void sparse_cholesky(compressed_matrix<T>* S, compressed_matrix<T>* L, compressed_matrix<T>* D);

#endif /* SPARSECHOLESKY_H_ */
