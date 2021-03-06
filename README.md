# matrixMult
Matrix operations are fundamental tools in Computer Graphics (CG). CG objects can be represented as sets of vertices (i.e., n-dimensional vectors with n being usually 4), and can be transformed (for example to animate them) through a composition of matrix operations. Matrix multiplication is the most common operation that is repeated millions of times when the object passes through the geometric pipeline. For efficiency, many of the operations are done concurrently. Nowadays, that is done by multi-core hardware on a video card, however, it all started as software libraries some years ago...

In this lab displays concurrent processing using threads by recreating the history: the task here was to implement multiplication of matrices concurrently computing the elements of the product.

Lets say that we want to multiply matrix A with the dimension m by k by a matrix B with the dimension k by n. To obtain the product, say matrix C, we compute dot products for each row of A and each column of B; i.e., to compute a specific element cij of the matrix C, we compute a dot product of row ai from the matrix A and the column bj from the matrix B.


The application should accept input in the following format:

the first line contains three numbers that specify the sizes of the matrices; for example, m, k, and n, followed by
m lines of k numbers that represent matrix A, and in turn followed by
k lines of n numbers that represent matrix B.
The program should assume the validity of m, k, and n (i.e., that the matrices A and B can be legally multiplied)
