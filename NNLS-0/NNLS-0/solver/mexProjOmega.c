/***********************************************************************
* ProjOmega: find the vector obtained from the matrix
*            U*V' by extracting the entries with indices specified
*            in (I,Jcol). 
* 
* xx = mexProjOmega(U,V,I,Jcol);
*  
* (I, Jcol) are in compressed sparse column format.
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

/********************************
* saxpymat:  z = z + alpha*y
* y dense matrix, z dense vector
********************************/
void saxpymat(const double *y,  const int idx1,
              const int istart, const int iend,
              const double alp, double *z)
{  int i;

   for(i=istart; i< iend-3; i++){             /* LEVEL 4 */
      z[i] += alp * y[i+idx1]; i++;
      z[i] += alp * y[i+idx1]; i++;
      z[i] += alp * y[i+idx1]; i++;
      z[i] += alp * y[i+idx1];
   }
   if(i < iend-1){                            /* LEVEL 2 */
      z[i] += alp * y[i+idx1]; i++;
      z[i] += alp * y[i+idx1]; i++;
   }
   if(i < iend){                              /* LEVEL 1 */
      z[i] += alp * y[i+idx1];
   }
   return;
}
/**********************************************************/
/**********************************************************/
void mexFunction(
      int nlhs,   mxArray  *plhs[], 
      int nrhs,   const mxArray  *prhs[] )
{    
     double   *U, *V, *irA, *jcA;
     double   *Vj, *xx, *Ut;
     int      isspU, isspV, m, n, nU, mV, len;   

     int      i, j, r, k, kstart, kend, rnU;
     double   tmp; 

/* CHECK THE DIMENSIONS */

   if (nrhs != 4) {
       mexErrMsgTxt("ProjOmega: requires 4 inputs"); }
   if (nlhs != 1) { 
       mexErrMsgTxt("ProjOmega: requires 1 output"); }

   /***** assign pointers *****/
       m  = mxGetM(prhs[0]); 
       nU = mxGetN(prhs[0]); 
       if (mxIsSparse(prhs[0])) { 
          mexErrMsgTxt("ProjOmega: U must be dense"); 
       } else {
          U = mxGetPr(prhs[0]); 
       }       
       mV = mxGetM(prhs[1]); 
       if (mxGetN(prhs[1]) != nU) { 
          mexErrMsgTxt("ProjOmega: U,V not compatible"); }     
       if (mxIsSparse(prhs[1])) { 
          mexErrMsgTxt("ProjOmega: V must be dense");      
       } else {
          V = mxGetPr(prhs[1]); 
       }       
       len = MAX(mxGetM(prhs[2]),mxGetN(prhs[2]));
       irA = mxGetPr(prhs[2]); 
       jcA = mxGetPr(prhs[3]); 
       n   = MAX(mxGetM(prhs[3]),mxGetN(prhs[3]))-1;
       /***** create return argument *****/
       plhs[0] = mxCreateDoubleMatrix(len,1,mxREAL); 
       xx = mxGetPr(plhs[0]); 

       /***** main body *****/
       Vj = mxCalloc(nU,sizeof(double));
       Ut = mxCalloc(m*nU,sizeof(double));
       for (j=0; j<m; j++) { 
	 for (k=0; k<nU; k++) {
	   Ut[k+j*nU] = U[j+k*m]; 
	 }
       }
       for (j=0; j<n; j++){
          for (k=0; k<nU; k++) { Vj[k] = V[j+k*mV]; }
          kstart = (int)jcA[j]; kend = (int)jcA[j+1]; 
          for (k=kstart; k<kend; k++) { 
	     r = (int)(irA[k]-1); /** adjust for matlab index **/ 
             rnU = r*nU; 
             tmp = 0; 
             for (i=0; i<nU; i++) { tmp += Vj[i]*Ut[i+rnU]; } 
             xx[k] = tmp; 
          }          
       }   
return;
}
/**********************************************************/

