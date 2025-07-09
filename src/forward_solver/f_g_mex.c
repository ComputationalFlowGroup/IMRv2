#include "mex.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    // Inputs: x, Ca, Ca1, Rst, l1, l2, v_a, v_nc, Rdot, Rstdot, mode (g or gdot)
    double *x = mxGetPr(prhs[0]);
    double Ca = mxGetScalar(prhs[1]);
    double Ca1 = mxGetScalar(prhs[2]);
    double Rst = mxGetScalar(prhs[3]);
    double l1 = mxGetScalar(prhs[4]);
    double l2 = mxGetScalar(prhs[5]);
    double v_a = mxGetScalar(prhs[6]);
    double v_nc = mxGetScalar(prhs[7]);
    double Rstdot = mxGetScalar(prhs[8]);
    int mode = (int) mxGetScalar(prhs[9]);

    size_t n = mxGetNumberOfElements(prhs[0]);

    // Outputs: g(x), gdot(x)
    plhs[0] = mxCreateDoubleMatrix(1, n, mxREAL);
   // plhs[1] = mxCreateDoubleMatrix(1, n, mxREAL);
   // double *gout = mxGetPr(plhs[0]);
   // double *gdotout = mxGetPr(plhs[1]);
    double *out = mxGetPr(plhs[0]);

    double Rst3m1 = pow(Rst, 3) - 1.0;
//    if (fabs(Rst3m1) < 1e-12)
//        Rst3m1 = (Rst3m1 < 0) ? -1e-12 : 1e-12;
    double Rst2 = Rst * Rst;

    for (size_t i = 0; i < n; ++i) {
        double xi = x[i];
        double xi2 = xi * xi;
        double xi3 = xi2 * xi;
        double xi3m1 = xi3 - 1.0;
        double xi5 = pow(xi,5);

        // Compute f_cy
        double common = pow(xi3m1/Rst3m1, 1.0/3.0);
        double num = l2 * common - 1.0;
        double den = 1.0 - l1 * common;
        double f = num / den;

        // Compute g(x)
        double term1 = 1.0/Ca;
        double term2a = (1.0/Ca1 - 1.0/Ca);
        double v_ncm1 = v_nc - 1;
        double gx = (1.0 / xi5 + 1.0 / xi2);
        if (mode == 1) {
        double term2 = term2a * pow(1 + pow(f, v_a), v_ncm1/v_a);
        double prefac = term1 + term2;
        // gout[i] = prefac * gx;
        out[i] = prefac * gx;
        }
        else if (mode == 2) {
        // Compute fdot_xy
        double fdnum = pow(xi3m1, 1.0/3.0);
        double fdden = pow(Rst3m1, 4.0/3.0) * den * den;
        double fd = fdnum * Rstdot * Rst2 * (l1 - l2) / fdden;
       
        // Compute gdot(x)
        double gdot = term2a * gx * v_ncm1 * pow(1 + pow(f, v_a), (v_ncm1 - v_a)/v_a) * pow(f,v_a-1) * fd;
        // gdotout[i] = gdot;
        out[i] = gdot;
        }
        else {
            mexErrMsgIdAndTxt("g_mex:invalidMode","Mode must be 1 (g) or 2 (gdot).");
        }
    }
}

