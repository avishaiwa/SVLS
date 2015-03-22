/***********************************************************************
* mexsvd.c : C mex file 
*
* Compilation in 
* unix: mex -O -largeArrayDims -lmwlapack -lmwblas mexsvd.c
*
* [U,S,VT,flag] = mexsvd(A,options); 
* A = U*S*VT 
*
* options = 1 (default), 
*             full SVD: 
*             dim(S) = M x N, dim(U) = M x M, dim(VT) = N x N, 
*         = 2, economical SVD: 
*              dim(S) = minMN x minMN, dim(U) = M x minMN, dim(VT) = minMN x N.
*
*  flag   = 0:  successful exit.
*         < 0:  if flag = -i, the i-th argument had an illegal value.
*         > 0:  LAPACK DBDSDC did not converge, updating process failed.
***********************************************************************/

#include <math.h>
#include <mex.h>
#include <matrix.h>
#include <string.h> /* needed for memcpy() */

#if !defined(MX_API_VER) || ( MX_API_VER < 0x07030000 )
typedef int mwIndex;
typedef int mwSize;
#endif

#if !defined(MAX)
#define  MAX(A, B)   ((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define  MIN(A, B)   ((A) < (B) ? (A) : (B))
#endif

/**********************************************************
* 
***********************************************************/
void mexFunction(
      int nlhs,   mxArray  *plhs[], 
      int nrhs,   const mxArray  *prhs[] )

{    double   *A, *U, *VT, *S, *flag, *work, *AA;  

     mwIndex  subs[2];
     mwSize   nsubs=2; 
     mwIndex  *irS, *jcS, *iwork; 
     mwSize   M, N, lwork, info, options, minMN, maxMN, k; 
     mwSize   LDA, LDU, LDVT, nU, mVT, mS, nS; 
     char     *jobz;

/* CHECK FOR PROPER NUMBER OF ARGUMENTS */

   if (nrhs > 2){
      mexErrMsgTxt("mexsvd: requires at most 2 input arguments."); }
   if (nlhs > 4){ 
      mexErrMsgTxt("mexsvd: requires at most 4 output argument."); }   

/* CHECK THE DIMENSIONS */

    M = mxGetM(prhs[0]); 
    N = mxGetN(prhs[0]); 
    if (mxIsSparse(prhs[0])) {
       mexErrMsgTxt("mexeig: sparse matrix not allowed."); }   
    A = mxGetPr(prhs[0]); 
    LDA = M; 
    minMN = MIN(M,N); 
    maxMN = MAX(M,N); 
    options = 1; 
    if (nrhs==2) { options = (int)*mxGetPr(prhs[1]); } 
    if (options==1) { 
       jobz="A"; nU = M; mVT = N; mS = M; nS = N; 
    } else {
       jobz="S"; nU = minMN; mVT = minMN; mS = minMN; nS = minMN; 
    }   
    LDU  = M; 
    LDVT = mVT; 
    /***** create return argument *****/
    plhs[0] = mxCreateDoubleMatrix(M,nU,mxREAL); 
    U = mxGetPr(plhs[0]);  
    plhs[1] = mxCreateSparse(mS,nS,minMN,mxREAL); 
    S   = mxGetPr(plhs[1]); 
    irS = mxGetIr(plhs[1]); 
    jcS = mxGetJc(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(mVT,N,mxREAL);     
    VT = mxGetPr(plhs[2]); 
    plhs[3] = mxCreateDoubleMatrix(1,1,mxREAL);     
    flag = mxGetPr(plhs[3]); 

    /***** Do the computations in a subroutine *****/   
    lwork = 4*minMN*minMN + MAX(maxMN,5*minMN*minMN+4*minMN);  
    work  = mxCalloc(lwork,sizeof(double)); 
    iwork = mxCalloc(8*minMN,sizeof(int)); 

    AA    = mxCalloc(M*N,sizeof(double)); 
    memcpy(AA,mxGetPr(prhs[0]),(M*N)*sizeof(double));

    dgesdd_(jobz,&M,&N, AA,&LDA,S,U,&LDU,VT,&LDVT,work,&lwork,iwork, &info); 

    for (k=0; k<minMN; k++) { irS[k] = k; }
    jcS[0] = 0;
    for (k=1; k<=nS; k++) { 
      if (k<minMN) { jcS[k] = k; } else { jcS[k] = minMN; }
    }  
    flag[0] = (float)info; 
    
    return;
 }
/**********************************************************/
