
[Introduction]
NNLS version 0: 
A MATLAB software for semidefinite programming
Copyright (c) 2009 by
Kim-Chuan Toh and Sangwoon Yun

This is software package for solving a nuclear norm
regularized least squares problem of the form: 

(P)  min_X {f(X) + mu*sum(svd(X))}
or 
(PS) min_X {f(X) + mu*trace(X) : X psd} 
where f(X) = 0.5*norm(A(X)-b)^2 + 0.5*rho*norm(X,'fro')^2.

The main algorithm is an acceleration proximal gradient 
method applied to (P); details can be found in the following 
reference: 

[Reference]
[1] K.C. Toh, and S.W. Yun    
    An accelerated proximal gradient algorithm for nuclear norm regularized 
    least squares problems,  preprint, National University of Singapore, 
    Apr 2009.  
Available at: http://www.optimization-online.org/DB_HTML/2009/03/2268.html

[Copyright] 
See Copyright.txt

--------------------------------------------------------------
[Installation] 
The software NNLS requires a few mex-files that you may need
to compile in MATLAB. 
To so, run MATALB in the directory NNLS, then type: 

>> Installmex 

After that, to see whether you have installed the software
correctly, type the following in MATLAB: 

>> runRandomMatComp
--------------------------------------------------------------