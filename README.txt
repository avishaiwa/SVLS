A MATLAB software for SVLS algorithm
Authors: Avishai Wagner and Or Zuk

This is software for solving problams of low matrix recovery
from rows and columns affine measurements using algorithm SVLS.

Installation:
=============

Download the SVLS directory from https://github.com/orzuk/SVLS. 
Use the "Download ZIP" button on the right hand side of the page. 
Once done, uncompress it, and follow the introductory files listed below.

You will have to change the path variable to your SVLS location. For example, if SVLS was saved under your home directory then enter the following in each introductory file: 
userDir = '~/SVLS/';


Dependencies: 
=============
The SVLS algorithm uses the PROPACKmod package for computing fast singular value decomposition. You can download PROPACKmod from here:  http://sun.stanford.edu/~rmunk/PROPACK/

Some of the files in the Simulations directory use other packages, since we compare SVLS to other algorithms. To be able to run these comparisons, you will have to download and install the following packages: 

1. optSpace, from : http://web.engr.illinois.edu/~swoh/software/optspace/code.html

2. SVT, from : http://svt.stanford.edu/code.html

3. APGL, from :http://www.math.nus.edu.sg/~mattohkc/NNLS.html

The Simulations directory also uses the package MaUtils for different basic usful functions.

Getting Started: 
================

An example running the SVLS algorithm is available in the file: test_SVLS

Directories: 
=============

Algorithms - recovery algorithms, including the SVLS algorithm

Simulations - Scripts for running and testing the package

Figs - directory with figures generating during the running of functions

Data - simulation results saved to files 


Files:
======
SVLS - implementation of the SVLS  algorithm 

Symmetric_SVLS - implement algorithm SVLS for symmetric matrices

gradient_descent - minimize the loss function (||Ar(X)-Br||_F)^2+(||(X)Ac-Bc||_F)^2

Symmetric_gradientDescent - minimize the loss function (||Ar(X)-Br||_F)^2

AffineMatrixRecovery - a 'Master' recovery function, recovers a matrix from affine measurements by calling the appropriate algorithm

samp_matrix - sample a random low-rank matrix 

All_Simulations - Generate random matrices and random measurements, estimate the original matrix using recovery algorithms and evaluate results 

make_SVLS_paper_figures - script running and producing all figures shown in paper [1].


Acknowledgment:
===============

SVLS was developed by Avishai Wagner and Or Zuk, as part of work on the papers,

[1] Low-Rank Matrix Recovery from Row-and-Column Affine Measurements", A. Wagner and O. Zuk

Please cite the above paper if using the package.

For support, any questions or comments, please contact:

Avishai Wagner: avishaiwa@gmail.com

Or Zuk: or.zuk@mail.huji.ac.il