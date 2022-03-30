///*
// * sparse_cholesky_matlab_api.cpp
// *
// *  Created on: Jul. 14, 2019
// *      Author: evangelia
// */

#include <chrono>
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
