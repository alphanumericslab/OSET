// `trimmed_filter_mex.cpp - MEX function wrapper for the trimmed_filter.cpp function (for compiling with MATLAB MEX command and running from MATLAB)

/*
% trimmed_filter - Trimmed filtering function.
%
%   [y] = trimmed_filter(x, mode, wlen, alpha, weights)
%
% Inputs:
%   x: Input signal (1D vector).
%   mode: String specifying the filter mode.
%         Available modes: 'mean', 'median', 'trmean', 'wmedian'.
%   wlen: Window length (integer).
%   alpha (optional): Alpha parameter (integer) used in 'trmean' mode.
%   weights (optional): Weights vector used in 'wmedian' mode.
%
% Outputs:
%   y: Filtered output (1D vector) of the same size as the input signal.
%
% Example:
%   x = [5.7, 7.5, 8, 8, 9, 17, -1, 22, 32, 1.5, 2, 9, -18, 7];
%   wlen = 5;
%   alpha = 2;
%   mode = 'trmean';
%   y = trimmed_filter(x, mode, wlen, alpha);
%
% References:
%   [1] Gonzalo R. Arce (2004). Nonlinear Signal Processing: A Statistical Approach.
%       John Wiley & Sons Inc.
%
% Revision History:
%   2008: First release
%   2023: Renamed from deprecated version TrimmedFilter.cpp
%
% Reza Sameni, 2022-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET
*/

#include <iostream>
#include <vector>
#include <cstring>
#include "list.h"
#include "trimmed_filter.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 3)
    {
        mexErrMsgTxt("Not enough input arguments. Usage: y = trimmed_filter(x, mode, wlen, alpha, weights)");
    }

    // Parse input data from MATLAB array
    double *x = mxGetPr(prhs[0]);
    int n = mxGetNumberOfElements(prhs[0]);
    char *type = mxArrayToString(prhs[1]);
    int w = (int)mxGetScalar(prhs[2]);
    int a = 0;
    double *h = nullptr;

    // Check for optional parameters
    if (nrhs >= 4)
    {
        a = (int)mxGetScalar(prhs[3]);
    }

    if (nrhs >= 5)
    {
        h = mxGetPr(prhs[4]);
    }

    // Check for valid input
    if (w < 1)
    {
        mexErrMsgTxt("Window length must be greater than 0.");
    }

    if (w > n)
    {
        mexErrMsgTxt("Window length exceeds the input vector length.");
    }

    // Check if the filter mode is valid
    std::string modeStr(type);
    if (modeStr != "mean" && modeStr != "median" && modeStr != "trmean" && modeStr != "wmedian")
    {
        mexErrMsgTxt("Invalid mode. Available modes: 'mean', 'median', 'trmean', 'wmedian'");
    }

    // Allocate memory for output array
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
    double *y = mxGetPr(plhs[0]);

    // Call trimmed_filter function
    int status = trimmed_filter(x, y, n, type, w, a, h);

    // Check for unknown filter type
    if (status)
    {
        mexErrMsgTxt("Unknown filter type.");
    }
}
