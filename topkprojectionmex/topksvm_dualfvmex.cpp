/*
 *  Description    : This file contains an MEX implementation to
 *                
 * Version History : 1.0, Sep. 13, 2017
 *
 * Author          : Dejun Chu, PhD Student, Tsinghua Unversity.
 *
 */

#include "mex.h"
#include <math.h>
#include <time.h>

#define dualfv_OUT plhs[0]

#define xtrn_IN prhs[0]
#define ytrn_IN prhs[1]
#define lambda_IN prhs[2]
#define A_IN prhs[3]
#define nY_IN prhs[4]
#define p_IN prhs[5]
#define num_IN prhs[6]

/*
 * compute the dual objective function value
 * xtrn: dim x num
 * ytrn: 1 x num
 * A   : nY x num, i.e., [alpha_1, alpha_2, ..., alpha_n]
 */

double topksvm_dualfv(double *xtrn, double *ytrn,
                      double lambda, double *A,
                      long nY, long p, long num)
{
    double dualfv = 0;

    double *Xiai, *sumXiai;
    Xiai = (double *)mxMalloc(sizeof(double)*(p*nY));
    sumXiai = (double *)mxMalloc(sizeof(double)*(p*nY));
    double *ptr_xi, *ptr_ai;
    double *ptr_Xiai = Xiai;
    double *ptr_sumXiai = sumXiai;

    double sumaici = 0, sumai;
    long i, j, k;
    long yi;
    double a_ij;  // A(j,i)
        
    for(i=0; i<p*nY; i++)
    {
        sumXiai[i] = 0;
        //mexPrintf("%f ", sumXiai[i]);
    }
    //mexPrintf("\n\n");

    for(i=0; i<num; i++)
    {
        yi = ytrn[i]-1;  // from matlab to C
        ptr_xi = xtrn + i*p;
        ptr_ai = A + i*nY;  

        sumai = 0;  // sum(alpha_i)
        for(j=0; j<nY; j++)
        {
            if(ptr_ai[j]!=0) // ==============>
            {
                a_ij = ptr_ai[j];
                sumai += a_ij;
            }
        }        

        sumaici = sumaici + sumai - ptr_ai[yi];  // 1st part of objective
        
        // mexPrintf("\ni: %d, ai'*ci: %f\n", i, sumai - ptr_ai[yi]);

        for(j=0; j<nY; j++)
        {
            if(ptr_ai[j]==0 && j!= yi) // ==============>
                continue;
            
            a_ij = ptr_ai[j];
            if(j==yi)
                a_ij = a_ij - sumai;

            ptr_Xiai = Xiai + j*p;
            ptr_sumXiai = sumXiai + j*p;
            
            //mexPrintf("\n i:%d, j:%d, Xiai:\n", i,j);            
            for(k=0; k<p; k++)
            {
                ptr_Xiai[k] = ptr_xi[k]*a_ij;
                ptr_sumXiai[k] += ptr_Xiai[k];
                
                //mexPrintf("    %f", ptr_Xiai[k]);
            }
        }
    }

    for(i=0; i<p*nY; i++)
        dualfv += sumXiai[i]*sumXiai[i]; 
    
    dualfv = -sumaici/num - dualfv/(2*num)/num/lambda;

    mxFree(sumXiai);
    mxFree(Xiai);

    return dualfv;
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{    
    long p, num;
    double *xtrn, *A, *dualfv;
    double *ytrn;
    double lambda;
    long nY;    
    
    if(nrhs!=7)
        mexErrMsgTxt("Wrong number of input arguments.");
    if(nlhs > 1 )
        mexErrMsgTxt("Too many output arguments.");
    
    xtrn = mxGetPr(xtrn_IN);
    ytrn = mxGetPr(ytrn_IN);
    A = mxGetPr(A_IN);
    lambda = mxGetScalar(lambda_IN);
    nY = mxGetScalar(nY_IN);
    p = mxGetScalar(p_IN);
    num = mxGetScalar(num_IN);
    
    /*mexPrintf("p: %d, num: %d\n", p, num);*/
        
    dualfv_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
    dualfv = mxGetPr(dualfv_OUT);
        
    dualfv[0] = topksvm_dualfv(xtrn, ytrn, lambda, A, nY, p, num);
    //mexPrintf("out_dualfv: %.12f\n", dualfv[0]);
    return;
}
