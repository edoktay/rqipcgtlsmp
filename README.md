# rqipcgtlsmp
MATLAB codes for performing RQI-PCGTLS-MP algorithm for solving total least squares problems.

This code can be used to reproduce the experiments in Mixed Precision Rayleigh Quotient Iteration for Total Least Squares Problems (http://arxiv.org/abs/2305.19028)


## Included MATLAB files
* **_pcgtls.m_** is a function that performs PCGTLS algorithm for RQI-PCGTLS and RQI-PCGTLS-MP with a user-specified precision.

* **_rqi.m_** is a function that performs RQI-PCGTLS (with uniform precision in fp64) and RQI-PCGTLS-MP (up to three precisions).

* **_alltests.m_** is an example script for comparing RQI-PCTLS and RQI-PCGTLS-MP (with three precisions) on random dense, \delta, Bjorck, Vanhuffel, and Toeplitz matrices defined in the manuscript.

* **_performance.m_** is a function that calculates speedup of RQI-PCGTLS-MP for various sized matrices and various numbers of Rayleigh quotient iterations performed in total.


## Requirements
* The codes have been developed and tested with MATLAB 2022a.
* The codes require some functions from libraries https://github.com/SrikaraPranesh/Multi_precision_NLA_kernels, https://github.com/higham/chop, and https://github.com/SrikaraPranesh/LowPrecision\_Simulation for half precision QR factorization.
* The codes require the Advanpix Multiprecision Computing Toolbox for extended precision computations. 
A free trial of Advanpix is available for download from https://www.advanpix.com/.


