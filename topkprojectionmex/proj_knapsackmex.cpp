/*
 *  File            : proj_knapsackmex.cpp
 *
 *  Version History : 0.1, Feb. 22, 2017
 *                    0.2, Oct. 9, 2017
 *
 *  Author          : Dejun Chu, PhD Student, Tsinghua Unversity.
 *
 *  Description     : This file contains an implementation in the C language
 *                    of algorithms described in the research paper:
 *
*/

#include "mex.h"
#include <math.h>
#include <time.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

#define X_OUT plhs[0]
#define t_OUT plhs[1]
#define A_IN prhs[0]
#define K_IN prhs[1]
#define R_IN prhs[2]
#define LAM0_IN prhs[3]
#define FIX_IN prhs[4]

#define POS_INF (-log(0.0))
#define NEG_INF (log(0.0))

/*
 *** compute funtion value and derivative of \varphi(lambda) ***
 */
double philambda(double *a, int dim,
                 double lambda, double r,
                 double li, double ui,
                 double *gleft, double *gright)
{
    double tmp, xi, pval = -r;
    double lb, ub;  // lower bound and upper bound
    double bdi = -1.0;  // bd = -b^2./d
    double deriv_left=0.0, deriv_right=0.0;

    for(int i=0; i<dim; i++)
    {
        tmp = a[i] - lambda;
        lb = a[i] - ui;
        ub = a[i] - li;

        if(tmp < li)
            xi = li;
        else if(tmp>=li && tmp<=ui)
            xi = tmp;
        else
            xi = ui;

        pval += xi;  // pval = \sum_i(bi*xi) - r

        if(lambda>lb && lambda<ub)
        {
            deriv_left += bdi;
            deriv_right += bdi;
        }
        if(lambda == lb)
            deriv_right += bdi;
        if(lambda == ub)
            deriv_left += bdi;
    }

    *gleft = deriv_left;
    *gright = deriv_right;

    return pval;
}

double philambda_fixing(double *a, double lambda, double r,
                        double li, double ui,
                        double *gleft, double *gright,
                        int *Iinx, int *Linx, int *Uinx,
                        int *num_Iinx, int *num_Linx, int *num_Uinx)
{
    double tmp, xi, pval;
    double lb, ub;  // lower bound and upper bound
    double bdi = -1.0;  // bd = -b^2./d
    double deriv_left=0.0, deriv_right=0.0;

    int inx, jI = 0, jL = 0, jU = 0;
    int numI = *num_Iinx;
    int numL = *num_Linx;  // num of all the breakpoints index in L set
    int numU = *num_Uinx;  // num of all the breakpoints index in U set

    pval = -r + numU*ui;  // bi = 1, pval = b'*x-r

    for(int i=0; i<numI; i++)
    {
        inx = Iinx[i];
        tmp = a[inx] - lambda;
        lb = a[inx] - ui;
        ub = a[inx] - li;

        if(tmp <= li)
        {
            xi = li;
            Linx[jL] = inx;
            jL++;
        }
        else if(tmp>li && tmp<ui)
        {
            xi = tmp;
            Iinx[jI] = inx;  // update Iinx set
            jI++;
        }
        else
        {
            xi = ui;
            Uinx[jU] = inx;
            jU++;
        }

        pval += xi;  // pval = \sum_i(bi*xi) - r

        if(lambda>lb && lambda<ub)
        {
            deriv_left += bdi;
            deriv_right += bdi;
        }
        if(lambda == lb)
            deriv_right += bdi;
        if(lambda == ub)
            deriv_left += bdi;
    }

    *gleft = deriv_left;
    *gright = deriv_right;

    if(pval>0)  // update Linx set
    {
        //for(int i=0; i<jL; i++)
        //    Linx[numL+i] = TLinx[i];
       numL = numL + jL;
       numI = numI - jL; // update the index number of Iinx set

       for(int i=0; i<jU; i++)
           Iinx[jI+i] = Uinx[i];
    }

    if(pval<0) // update Uinx set
    {
        //for(int i=0; i<jU; i++)
        //    Uinx[numU+i] = TUinx[i];
        numU = numU + jU;
        numI = numI - jU; // update the index number of Iinx set

        for(int i=0; i<jL; i++)
            Iinx[jI+i] = Linx[i];
    }

    *num_Iinx = numI;
    *num_Linx = numL;
    *num_Uinx = numU;

    return pval;
}


/*
 *** find the smalllest breakpoint in (alpha, beta) ***
 */
int findmin(double *bps, int *num,
            double alpha, double beta,
            double *min)
{
    int j = 0;
    double minimum = POS_INF;

    for(int i=0; i< *num; i++)
    {
        if(bps[i]>alpha && bps[i]<beta)
        {
            bps[j] = bps[i];
            j++;

            if(bps[i]<minimum)
                minimum = bps[i];
        }
    }

    *num = j;
    *min = minimum;

    if(minimum == POS_INF)
        return 1;
    else
        return 0;

}

/*
 *** find the largest breakpoint in (alpha, beta) ***
 */
int findmax(double *bps, int *num,
            double alpha, double beta,
            double *max)
{
    int j = 0;
    double maximum = NEG_INF;

    for(int i=0; i< *num; i++)
    {
        if(bps[i]>alpha && bps[i]<beta)
        {
            bps[j] = bps[i];
            j++;

            if(bps[i]>maximum)
                maximum = bps[i];
        }
    }

    *num = j;
    *max = maximum;

    if(maximum == NEG_INF)
        return 1;
    else
        return 0;

}

/* Solve a knapsack problem
 *   phi(lambda) =
 *      sum_i{ bi MIN{MAX(li, (-lam bi+ai)/di) ,ui} } - r
 *
 *when s=r, f(t)=
 *      sum_i{ MIN{MAX(0, ai-t), r/k} } - r
 *
 * So
 *   bi = 1, ai = ai, di = 1, lam = t;, li = 0, ui = r/k.
 * */

double newtonknapsack(double *a, double *x,
                   int dim, int k,
                   double r, double tol,
                   double lam0, int FIX)
{
    double bi, di, li, ui;
    double lambda_k, lambda_N;
    double alpha_k, alpha_old, beta_k, beta_old;
    alpha_old = NEG_INF;
    beta_old = POS_INF;

    bi = 1.0;
    di = 1.0;
    li = 0;
    ui = r/k;

    double *bps; // bps: breakpoints in (alpha, beta)
    double bps_min = POS_INF, bps_max = NEG_INF;
    int num_bps = 2*dim; // number of breakpoints in (alpha, beta)
    bps = (double*)mxMalloc(num_bps*sizeof(double));

    int *Iinx, *Linx, *Uinx; // save index for fixing variable
    int num_Iinx = dim, num_Linx = 0, num_Uinx = 0; // index num of I, L and U sets

    if(FIX)
    {
        Iinx = (int*)mxMalloc(dim*sizeof(int));
        Linx = (int*)mxMalloc(dim*sizeof(int));  // save temporary index from Iinx set
        Uinx = (int*)mxMalloc(dim*sizeof(int));  // that fixed in the iteration

        for(int i=0; i<dim; i++)
            Iinx[i] = i;  // initializing the index
    }

    double tmp;
    int iter;

    /* phi_lamk is the function value of \varphi(\lambda_k)
     * gleft is the left derivative of \varphi'(\lambda_k)
     * */
    double phi_lamk, gleft = 0.0, gright = 0.0;
    double phi_alphak = 0.0, phi_betak = 0.0;
    double phi_alphaold, phi_betaold;
    phi_alphaold = dim*ui - r;
    phi_betaold = -r;  // dim*li-r

    double slope, lambda_S; //used in secant step

    /* set tolerence to avoid the huge lambda
     * when derivative is a small negative
     */
    double grad_tolerence = -0.001;

    iter = 0;
    lambda_k = lam0;


    for(int i=0; i<dim; i++)
    {
        tmp = a[i]-ui;  // (a[i]-di*ui)/bi
        bps[i] = tmp;
        // find the smallest breakpoint in (alpha_k, lambda_k)
        if(tmp<lambda_k && tmp<bps_min)
            bps_min = tmp;

        tmp = a[i];     // (a[i]-di*li)/bi
        bps[dim+i] = tmp;
        // find the largest breakpoint in (lambda_k, beta_k)
        if(tmp>lambda_k && tmp>bps_max)
            bps_max = tmp;
    }

    phi_lamk = philambda(a, dim, lambda_k, r, li, ui, &gleft, &gright);

    while(fabs(phi_lamk)>tol)
    {
        iter++;

        if(phi_lamk>0)
        {
            alpha_k = lambda_k;
            alpha_old = alpha_k;
            beta_k = beta_old;

            phi_alphak = phi_lamk;
            phi_alphaold = phi_alphak;
            phi_betak = phi_betaold;

            if(gright < grad_tolerence)  // grad_right<0
            {
                lambda_N = lambda_k - phi_lamk/gright;
                if(lambda_N < beta_k)
                    lambda_k = lambda_N;
                else // lambda_N jumps out of [alpha_k, beta_k]
                {
                    slope = (phi_betak-phi_lamk)/(beta_k - lambda_k);
                    lambda_S = lambda_k - phi_lamk/slope;
                    lambda_k = lambda_S<bps_max ? lambda_S:bps_max;
                }
            }
            else  // grad_right==0
            {   // find the smallest breakpoint in (alpha_k, beta_k)
                int flag;
                flag = findmin(bps, &num_bps,
                               alpha_k, beta_k, &bps_min);
                if(flag)
                {
                    printf("Infeasible problem!\n");
                    return flag;  //Will the projection is not feasible？
                }
                else
                    lambda_k = bps_min;
            }
        }
        else // phi_lamk<0
        {
            alpha_k = alpha_old;
            beta_k = lambda_k;
            beta_old = beta_k;

            phi_alphak = phi_alphaold;
            phi_betak = phi_lamk;
            phi_betaold = phi_betak;

            if(gleft < grad_tolerence)  //grad_left<0
            {
                lambda_N = lambda_k - phi_lamk/gleft;
                if(lambda_N > alpha_k)
                    lambda_k = lambda_N;
                else  // secant step
                {
                    slope = (phi_alphak - phi_lamk)/(alpha_k - lambda_k);
                    lambda_S = lambda_k - phi_lamk/slope;
                    lambda_k = lambda_S>bps_min ? lambda_S:bps_min;
                }
            }
            else  // grad_left == 0
            {   //find the largest breakpoint in (alpha_k, beta_k)
                int flag;
                flag = findmax(bps, &num_bps,
                               alpha_k, beta_k, &bps_max);
                if(flag)
                {
                    printf("Infeasible problem!\n");
                    return flag;  //Will the projection is not feasible？
                }
                else
                    lambda_k = bps_max;
            }
        }

        if(FIX)  // using variable fixing stragety
            phi_lamk = philambda_fixing(a, lambda_k, r,
                                        li, ui, &gleft, &gright,
                                        Iinx, Linx, Uinx,
                                        &num_Iinx, &num_Linx, &num_Uinx);
        else
            phi_lamk = philambda(a, dim, lambda_k, r, li, ui, &gleft, &gright);
    }

    for(int i=0; i<dim; i++)
    {
        tmp = a[i] - lambda_k;
        if(tmp<li)
            x[i] = li;
        else if(tmp>ui)
            x[i] = ui;
        else
            x[i] = tmp;
    }

    if(FIX)
    {
        mxFree(Uinx);
        mxFree(Linx);
        mxFree(Iinx);
    }

    mxFree(bps);

    return lambda_k;
}

// [a_topk, t] = proj_knapsackmex(a, k, r, t0, FIX);
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double *a, *x, *t;
    int k, dim;
    double r, lam0, tol = 1e-8;
    int FIX;

    if(nrhs<1 || nrhs>5)
        mexErrMsgTxt("Wrong number of input arguments.\n Usage:"
                     "[a_topk, t] = proj_knapsackmex(a, k, r, t0, FIX).\n");
    if(nlhs > 2 )
        mexErrMsgTxt("Too many output arguments.\n Usage:"
                     "[a_topk, t] = proj_knapsackmex(a, k, r, t0, FIX).\n");

    switch(nrhs)
    {
    case 1:
        k = 1;
        r = 1.0;
        lam0 = 0;
        FIX = 0;
        break;
    case 2:
        k = mxGetScalar(K_IN);
        r = 1.0;
        lam0 = 0;
        FIX = 0;
        break;
    case 3:
        k = mxGetScalar(K_IN);
        r = mxGetScalar(R_IN);
        lam0 = 0;
        FIX = 0;
        break;
    case 4:
        k = mxGetScalar(K_IN);
        r = mxGetScalar(R_IN);
        lam0 = mxGetScalar(LAM0_IN);
        FIX = 0;
        break;
    case 5:
        k = mxGetScalar(K_IN);
        r = mxGetScalar(R_IN);
        lam0 = mxGetScalar(LAM0_IN);
        FIX = mxGetScalar(FIX_IN);
        break;
    default:
        mexErrMsgTxt("Wrong number of input arguments.");
    }
    //mexPrintf("Input: k=%d, r=%f, lam0=%f, FIX=%d\n", k, r, lam0, FIX);
    
    dim = mxGetM(A_IN);
    a = mxGetPr(A_IN);
    X_OUT = mxCreateDoubleMatrix(dim, 1, mxREAL);
    x = mxGetPr(X_OUT);
    t_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
    t = mxGetPr(t_OUT);
        
    *t = newtonknapsack(a, x, dim, k, r, tol, lam0, FIX);
    return;
}
