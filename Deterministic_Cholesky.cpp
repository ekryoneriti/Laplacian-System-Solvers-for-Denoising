/*
 * Deterministic_Cholesky.cpp
 *
 *  Created on: Jul. 11, 2019
 *      Author: evangelia
 */
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/storage.hpp>



//int main () {
//    using namespace boost::numeric::ublas;
//    matrix<double> m (3, 3);
//    for (unsigned j = 0; j < m.size2 (); ++ j) {
//        matrix_column<matrix<double> > mc (m, j);
//        for (unsigned i = 0; i < mc.size (); ++ i)
//            mc (i) = 3 * i + j;
//        std::cout << mc << std::endl;
//    }
//}
