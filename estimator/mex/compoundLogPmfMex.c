/* compoundLogPmfMex.c  --  log P(Z = z), z = 0..zmax, for a compound Poisson Z.
 *
 *   lp = compoundLogPmfMex(lambda, jumps, zmax)
 *
 *   lambda : intensities (vector, >= 0)
 *   jumps  : jump sizes  (vector of positive integers, same length as lambda)
 *   zmax   : largest count (nonnegative integer scalar)
 *   lp     : (zmax+1) x 1 column, lp(z+1) = log P(Z = z)
 *
 * Panjer recursion in the log domain, term by term as in the MATLAB fallback of
 * estimateMuEfficientGMM>compoundLogPmf (same coefficients, same order of
 * summation), so the two agree to rounding:
 *   lp(0) = -sum(lambda)
 *   lp(z) = logsumexp_k [ log(k*lambda_k) + lp(z-k) ] - log(z)
 *
 * Build with buildPaperMex (writes the binary outside the library).
 */
#include <math.h>
#include <stdlib.h>

#ifndef COMPOUND_NO_MEX
#include "mex.h"
#endif

/* Core recursion; lp must hold zmax+1 values. Returns 0 on success. */
int compoundLogPmfCore(const double *lambda, const double *jumps, long K,
                       long zmax, double *lp)
{
    long k, z, kmax = 0;
    double tot = 0.0;
    double *coef, *lc;
    int *pos;
    for (k = 0; k < K; k++) {
        long j = (long) jumps[k];
        if (j > kmax) kmax = j;
        tot += lambda[k];
    }
    for (z = 0; z <= zmax; z++) lp[z] = -INFINITY;
    lp[0] = -tot;
    if (kmax < 1 || zmax < 1) return 0;
    coef = (double *) calloc((size_t) kmax + 1, sizeof(double));
    lc   = (double *) calloc((size_t) kmax + 1, sizeof(double));
    pos  = (int *)    calloc((size_t) kmax + 1, sizeof(int));
    if (!coef || !lc || !pos) { free(coef); free(lc); free(pos); return 1; }
    /* coef(jumps) = jumps .* lambda  (as in MATLAB: a later duplicate overwrites) */
    for (k = 0; k < K; k++) {
        long j = (long) jumps[k];
        if (j >= 1) coef[j] = (double) j * lambda[k];
    }
    for (k = 1; k <= kmax; k++) {
        pos[k] = coef[k] > 0.0;
        lc[k] = pos[k] ? log(coef[k]) : -INFINITY;
    }
    for (z = 1; z <= zmax; z++) {
        long top = kmax < z ? kmax : z;
        double m = -INFINITY, acc = 0.0, t;
        int any = 0;
        for (k = 1; k <= top; k++) {                 /* largest term */
            if (!pos[k]) continue;
            t = lc[k] + lp[z - k];
            if (t > -INFINITY) {
                if (!any || t > m) m = t;
                any = 1;
            }
        }
        if (!any) continue;                          /* lp[z] stays -Inf */
        for (k = 1; k <= top; k++) {                 /* sum of exp(term - largest) */
            if (!pos[k]) continue;
            t = lc[k] + lp[z - k];
            if (t > -INFINITY) acc += exp(t - m);
        }
        lp[z] = m + log(acc) - log((double) z);
    }
    free(coef); free(lc); free(pos);
    return 0;
}

#ifndef COMPOUND_NO_MEX
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    long K, zmax;
    double zm;
    (void) nlhs;
    if (nrhs != 3)
        mexErrMsgIdAndTxt("compoundLogPmfMex:Inputs", "Usage: lp = compoundLogPmfMex(lambda,jumps,zmax)");
    if (!mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) ||
        mxIsComplex(prhs[0]) || mxIsComplex(prhs[1]))
        mexErrMsgIdAndTxt("compoundLogPmfMex:Type", "Inputs must be real double.");
    K = (long) mxGetNumberOfElements(prhs[0]);
    if ((long) mxGetNumberOfElements(prhs[1]) != K)
        mexErrMsgIdAndTxt("compoundLogPmfMex:Size", "lambda and jumps must have the same length.");
    zm = mxGetScalar(prhs[2]);
    if (!(zm >= 0) || zm != floor(zm))
        mexErrMsgIdAndTxt("compoundLogPmfMex:Zmax", "zmax must be a nonnegative integer.");
    zmax = (long) zm;
    plhs[0] = mxCreateDoubleMatrix((mwSize) zmax + 1, 1, mxREAL);
    if (compoundLogPmfCore(mxGetPr(prhs[0]), mxGetPr(prhs[1]), K, zmax, mxGetPr(plhs[0])) != 0)
        mexErrMsgIdAndTxt("compoundLogPmfMex:Memory", "Out of memory.");
}
#endif
