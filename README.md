The scrips in this repository are numerical implementations of the paper (tomp_siam.pdf)

**Thresholding Orthogonal Matching Pursuit for Imaging Applications**, *Hai Le, Alexei Novikov*.

The *main.m* file (the only file to run) compares three algorithms: SquareRoot_LASSO, Stagewise Orthogonal Matching Pursuit (StOMP), and Thresholding Orthogonal Matching Pursuit (TOMP) for the sparse recovery problem:

$Ax=y$

where A is a measurement matrix, y is a data vector, and x is the sparse vector that we want to recover.

The settings for the Ax = y come from the settings of the paper https://arxiv.org/abs/1908.01479 where A is a highly coherent matrix (nearby normalized column vectors are almost parallel to each other).
