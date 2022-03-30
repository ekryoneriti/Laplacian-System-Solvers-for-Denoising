/*
 * SparceCholesky.h
 *
 *  Created on: Jul. 14, 2019
 *      Author: evangelia
 */



#ifndef SPARSECHOLESKY_MAT_H_
#define SPARSECHOLESKY_MAT_H_

#include <boost/numeric/ublas/matrix_sparse.hpp>

using namespace boost::numeric::ublas;

template <typename T>
void sparse_cholesky(compressed_matrix<T>* S, compressed_matrix<T>* L, compressed_matrix<T>* D);

#endif /* SPARSECHOLESKY_MAT_H_ */
