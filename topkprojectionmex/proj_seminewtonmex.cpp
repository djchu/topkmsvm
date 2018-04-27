/*
 *  Description     : This file contains an implementation in the C language
 *                    of algorithms described in the research paper:
 *
 *                    It solves the biased top-k simplex projection problem
 *                    with semi-smooth newton method.
 *
 * File            : proj_seminewtonmex.cpp
 * Version History : 1.0, Jan. 12, 2017
 * Version History : 1.1, Mar. 27, 2017
 *
 * Author          : Dejun Chu, PhD Student, Tsinghua Unversity.
 *
 */

#include "mex.h"
#include <math.h>
#include <time.h>

#define X_OUT plhs[0]
#define T_OUT plhs[1]
#define J_OUT plhs[2]  // save the index j if a[j]>t_max
#define Jlength_OUT plhs[3]

#define A_IN prhs[0]
#define K_IN prhs[1]
#define RHO_IN prhs[2]

double sumklargest(double *a, int dim, int k, double *a_kth)
{
    int i, k_inx, iter = 0, Ulength, gi, li;
    double sumk = 0.0, a_k;
    unsigned int *U, *G, *L;
    U = (unsigned int *)mxMalloc(dim*sizeof(unsigned int));
    G = (unsigned int *)mxMalloc(dim*sizeof(unsigned int));
    L = (unsigned int *)mxMalloc(dim*sizeof(unsigned int));
    
    Ulength = dim;
    for(int j=0; j<dim; j++)
        U[j] = j;
    
    /*mexPrintf("k_inx: ");*/
    srand(time(NULL));
    while(1)
    {
        iter++;
        
        i = rand() % Ulength; // i \in [0,...|U|-1]
        k_inx = U[i];
        a_k = a[k_inx];
        
        /*mexPrintf("%d ", k_inx);*/
        
        gi = 0;  /* gi is the number of G */
        li = 0;  /* li is the number of L */
        for(int j=0; j<Ulength; j++)
        {
            if(a[U[j]] < a_k)
            {
                L[li] = U[j];
                li++;
            }else
            {
                G[gi] = U[j];
                gi++;
            }
        }
        
        if(gi==k+1)
        {
            double minG = 1.0/0.0;
            
            for(int j=0; j<gi; j++)
            {
                sumk += a[G[j]];
                if(a[G[j]]<=minG & G[j]!=k_inx)
                    minG = a[G[j]];
            }
            
            sumk -= a_k;
            *a_kth = minG;
            return sumk;
        }
        
        if(gi>=k)
        {
            if(gi==k)
            {
                for(int j=0; j<gi; j++)
                    sumk += a[G[j]];
                
                *a_kth = a_k;
                return sumk;
            }
            else
            {
                int kk = 0;
                for(int j=0; j<gi; j++)
                {
                    if(G[j]==k_inx)
                        continue;
                    U[kk] = G[j];
                    kk++;
                }
                Ulength = gi-1;
            }
        }
        else
        {
            for(int j=0; j<gi; j++)
                sumk += a[G[j]];
            
            k = k-gi;
            
            for(int j=0; j<li; j++)
                U[j] = L[j];
            
            Ulength = li;
        }
    }
    /*mexPrintf("\n");*/
    
    mxFree(L);
    mxFree(G);
    mxFree(U);
    return sumk;
}


void jacobi(double *a, int dim, double s, double t,
            double k, double rho,
            double *F_ts,
            double *detJ, double *invJ)
{
    int i, numI1 = 0, numI2 = 0;
    double ai;
    double fval, gval;
    double ftsgrad_t, ftsgrad_s;
    double gstgrad_t, gstgrad_s;
    double eps = 1.0e-12;
    
    fval = -s;
    gval = k*t - k*rho*s;
    
    for(i=0; i<dim; i++)
    {
        ai = a[i];
        double ai_t_sk = ai - t - s/k;
        double ai_t = ai - t;
        
        if(ai_t_sk>eps)  /* a[i]-t>s/k */
        {
            numI1 += 1;
            fval += s/k;
            gval += ai_t_sk;
            
        }else
            if(ai_t>0)  /* 0<a[i]-t<=s/k */
            {
                numI2 += 1;
                fval += ai_t;
            }
    }
    
    F_ts[0] = fval;
    F_ts[1] = gval;
    
    ftsgrad_t = -numI2;
    ftsgrad_s = numI1/k - 1;
    gstgrad_t = -numI1 + k;
    gstgrad_s = -numI1/k - k*rho;
    
    *detJ = ftsgrad_t*gstgrad_s - gstgrad_t*ftsgrad_s;
    if(*detJ == 0)
        return;
    else
    {
        invJ[0] = gstgrad_s / (*detJ);
        invJ[1] = -ftsgrad_s / (*detJ);
        invJ[2] = -gstgrad_t / (*detJ);
        invJ[3] = ftsgrad_t / (*detJ);
    }
}

void init_ts(double *a,
        int dim,
        int k,
        double rho,
        double *t,
        double *s)
{
    double sumk;
    double a_kth;
    sumk = sumklargest(a, dim, k, &a_kth);
    
    *s = sumk/(1+k*rho);
    *t = a_kth - *s/k;
}

double seminewton(double *a, double *x, int dim,
                  int k, double rho, double tol)
{
    int iter = 0;
    double s0;
    double t = 0.11, s = 0.22;
    
    double detJ;
    double *F, *invJ, *d;
    
    F = (double*)mxMalloc(2*sizeof(double));
    invJ = (double*)mxMalloc(4*sizeof(double));
    d = (double*)mxMalloc(2*sizeof(double));
    
    init_ts(a, dim, k, rho, &t, &s);
    s0 = s;
    
    /*
     * mexPrintf("Newton-Raphson algorithm begins...\n");
     * mexPrintf(" iter: %d, t: %f, s: %f\n",
     * iter, t, s);
     * fflush(stdout);*/
    
    if(s0<=0)
    {
        s = 0;
        t = 0;
    }
    else
    {
        double error = 2*tol;
        while(error > tol)
        {
            iter = iter+1;
            jacobi(a, dim, s, t, k, rho, F, &detJ, invJ);
            
            if(detJ == 0)
            {
                s = s0;
                break;
            }
            else
            {
                d[0] = invJ[0]*F[0] + invJ[1]*F[1];
                d[1] = invJ[2]*F[0] + invJ[3]*F[1];
                t -= d[0];
                s -= d[1];
            }
            
            error = F[0]*F[0] + F[1]*F[1];
            
            /*
            mexPrintf(" iter: %d, t: %f, s: %f, error:%f\n",
               iter, t, s, error);
            fflush(stdout); */
        }
    }
    /*
     * mexPrintf("Newton-Raphson algorithm has been finished.\n");
     * mexPrintf("|J|: %f, t: %f, s: %f.\n", detJ, t, s);
     * fflush(stdout); */
    
    for(int i=0; i<dim; i++)
    {
        double ai_t_sk = a[i]-t-s/k;
        double ai_t = a[i] - t;
        x[i] = 0;
        
        if(ai_t_sk>0)
            x[i] = s/k;
        else
            if(ai_t>0)
                x[i] = ai_t;
    }
    
    mxFree(d);
    mxFree(invJ);
    mxFree(F);
    
    return t;
}


/*
 * Reduced version of seminewton where the element a_j is
 * deleted from a if a_j<=t_max
 */
double seminewton_reduced(double *a, double *x, int dim,
                          double *j_gtmax, double *jlength,
                          int k, double rho, double tol)
{
    int iter = 0;
    double s0;
    double t = 0.11, s = 0.22;

    double detJ;
    double *F, *invJ, *d;

    F = (double*)mxMalloc(2*sizeof(double));
    invJ = (double*)mxMalloc(4*sizeof(double));
    d = (double*)mxMalloc(2*sizeof(double));

    init_ts(a, dim, k, rho, &t, &s);
    s0 = s;

    /*
     * mexPrintf("Newton-Raphson algorithm begins...\n");
     * mexPrintf(" iter: %d, t: %f, s: %f\n",
     * iter, t, s);
     * fflush(stdout);*/

    if(s0<=0)
    {
        s = 0;
        t = 0;
    }
    else
    {
        // ******* Delete a_j from a  **********
        int i = 0;
        for(int j=0; j<dim; j++)
        {
            if(a[j]>t) // a_j <= t_max
            {
                if(j>i)
                    a[i] = a[j];

                j_gtmax[i] = j+1;  // matlab has 1-based indexing
                i++;
            }
        }
        dim = i;
        *jlength = dim;
        // ******* Delete a_j from a  **********

        double error = 2*tol;
        while(error > tol)
        {
            iter = iter+1;
            jacobi(a, dim, s, t, k, rho, F, &detJ, invJ);

            if(detJ == 0)
            {
                s = s0;
                break;
            }
            else
            {
                d[0] = invJ[0]*F[0] + invJ[1]*F[1];
                d[1] = invJ[2]*F[0] + invJ[3]*F[1];
                t -= d[0];
                s -= d[1];
            }

            error = F[0]*F[0] + F[1]*F[1];

            /*
            mexPrintf(" iter: %d, t: %f, s: %f, error:%f\n",
               iter, t, s, error);
            fflush(stdout); */
        }
    }
    /*
     * mexPrintf("Newton-Raphson algorithm ends.\n");
     * mexPrintf("|J|: %f, t: %f, s: %f.\n", detJ, t, s);
     * fflush(stdout); */

    for(int i=0; i<dim; i++)
    {
        double ai_t_sk = a[i]-t-s/k;
        double ai_t = a[i] - t;
        x[i] = 0;

        if(ai_t_sk>0)
            x[i] = s/k;
        else
            if(ai_t>0)
                x[i] = ai_t;
    }

    mxFree(d);
    mxFree(invJ);
    mxFree(F);

    return t;
}


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double *a, *x;
    double *j_gtmax, *jlength;
    int k, dim;
    double rho, tol = 1e-6;
    double *t_newton;
    
    if(nrhs<1 || nrhs>3)
        mexErrMsgTxt("Wrong number of input arguments.\nUsage:\n"
                     "  [z_topk, t] = proj_seminewtonmex(a, k ,r).\n"
                     "  or \n"
                     "  [z_topk, t, j_gt, jlength] = proj_seminewtonmex(a, k ,r).\n");
    if(nlhs > 4 )
        mexErrMsgTxt("Too many output arguments.\nUsage:\n"
                     "  [z_topk, t] = proj_seminewtonmex(a, k ,r).\n"
                     "  or \n"
                     "  [z_topk, t, j_gt, jlength] = proj_seminewtonmex(a, k ,r).\n");

    k = mxGetScalar(K_IN);
    rho = mxGetScalar(RHO_IN);
    
    dim = mxGetM(A_IN);
    a = mxGetPr(A_IN);
    X_OUT = mxCreateDoubleMatrix(dim, 1, mxREAL);
    x = mxGetPr(X_OUT);
    T_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
    t_newton = mxGetPr(T_OUT);
    
    if(nlhs<3)
        *t_newton = seminewton(a, x, dim, k, rho, tol);
    else
    {
        //mexPrintf("Reduced seminewton method.\n");
        J_OUT = mxCreateDoubleMatrix(dim, 1, mxREAL);
        j_gtmax = mxGetPr(J_OUT);
        Jlength_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
        jlength = mxGetPr(Jlength_OUT);
        *jlength = 0;

        *t_newton = seminewton_reduced(a, x, dim, j_gtmax, jlength,
                                       k, rho, tol);
    }
    return;
}
