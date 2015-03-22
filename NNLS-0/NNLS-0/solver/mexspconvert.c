/***********************************************************************
* mexspconvert: form a sparse matrix from the V with indices 
*               recorded in (I,Jcol). 
* 
* xx = mexspconvert(nr,nc,V,I,Jcol);
*  
* (I,Jcol) are in compressed sparse column format.
*
* NNLS, version 0: 
* Copyright (c) 2009 by
* Kim-Chuan Toh and Sangwoon Yun 
************************************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>

#if !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 )
typedef int mwIndex;
typedef int mwSize;
#endif

#if !defined(MAX)
#define  MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

/**********************************************************/
void mexFunction(
      int nlhs,   mxArray  *plhs[], 
      int nrhs,   const mxArray  *prhs[] )
{    
     double   *A, *V, *I, *Jcol; 
     mwIndex  *irA, *jcA;
     int      m, n, nnz, j, k;   

/* CHECK THE DIMENSIONS */

       if (nrhs != 5) {
          mexErrMsgTxt("mexspconvert: requires 5 inputs"); }
       if (nlhs != 1) { 
          mexErrMsgTxt("mexspconvert: requires 1 output"); }
/***** assign pointers *****/
       m  = (int)*mxGetPr(prhs[0]); 
       n  = (int)*mxGetPr(prhs[1]); 
       if (mxIsSparse(prhs[2])) { 
          mexErrMsgTxt("mexspconvert: V must be dense");      
       } else {
          V = mxGetPr(prhs[2]); 
       }   
       nnz = mxGetM(prhs[2]);
       if (mxGetM(prhs[3]) != nnz) {
          mexErrMsgTxt("mexspconvert: V and I not compatible");
       }
       I =  mxGetPr(prhs[3]); 
       if (MAX(mxGetM(prhs[4]),mxGetN(prhs[4])) != n+1) {
          mexErrMsgTxt("mexspconvert: Jcol not compatible");
       } 
       Jcol = mxGetPr(prhs[4]); 
       /***** create return argument *****/
       plhs[0] = mxCreateSparse(m,n,nnz,mxREAL); 
       A   = mxGetPr(plhs[0]); 
       irA = mxGetIr(plhs[0]); 
       jcA = mxGetJc(plhs[0]); 

       /***** main body *****/
       jcA[0] = 0;  
       for (k=0; k<nnz; k++){         
           A[k]   = V[k];
           irA[k] = (int)(I[k]-1); /** adjust for matlab index **/ 
       } 
       for (j=0; j<n; j++){         
          jcA[j+1] = Jcol[j+1];
       }
return;
}
/**********************************************************/

