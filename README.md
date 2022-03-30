# Laplacian-System-Solvers-for-Denoising

Recently there has been a breakthrough in solving Laplacian linear systems [2, 3]. Laplacian linear systems are related with the applications of Laplacian 
denoising and Total Variation denoising. This is because these problems reduce to solving (a sequence of) banded Laplacian linear systems. The new 
Laplacian system solvers have asymp- totic running time O(n2 log n2), where n2 is the number of pixels of the given image, which is better than standard 
LU factorization with asymptotic running time O(n3). To the best of my knowledge, the work on Laplacian linear system solvers is theoretical and there is 
no work on empirical performance. My goal in this project is to compare the empirical running time of the new Laplacian system solvers against standard 
numerical linear algebra solvers such as LU (with pivoting) and Conjugate Gradients for the applications of Laplacian denoising and Total Variation 
denoising. For the former application, the problem reduces to a single banded Laplacian linear system but for the latter application the problem reduces 
to non-linear equations which we will solve using Newton’s method. At each iteration of Newton’s method we will have to solve a banded Laplacian linear 
system. 

The published papers that I chose to use for the project are “Approximate Gaussian Elimination for Laplacians - Fast, Sparse and Simple” [1] and the 
review paper “Laplacian Solvers and Their Algorithmic Applications” [4]. The latter paper describes a sparse approximate Gaussian elimination algorithm 
for Laplacian matrices. This is the algorithm with best known asymptotic running time for solving Laplacian linear systems.

References
[1] Rasmus Kyng and Sushant Sachdeva. Approximate Gaussian elimination for Laplacians - Fast, Sparse, and Simple. 2016 IEEE 57th Annual Symposium on 
Foundations of Computer Science (FOCS), 2016.
[2] D. A. Spielman and S.-H. Teng. Nearly-linear time algorithms for graph partition- ing, graph sparsification, and solving linear systems. ACM Symposium 
on Theory of Computing (STOC), 2004.
[3] D. A. Spielman and S.-H. Teng. Nearly-linear time algorithms for preconditioning and solving symmetric, diagonally dominant linear systems. 
SIAM J. MATRIX ANAL. APPL, 35(3):835–885, 2006.
[4] Nisheeth K. Vishnoi. Lx b Laplacian Solvers and Their Algorithmic Applications. Foundations and Trends in Theoretical Computer Science, 2013.
