// resampler_mex.cpp: MATLAB mexFunction interface for the Resampler class

#include <iostream>
#include "mex.h"
#include "resampler.h"

// Function to perform resampling using the Resampler class
void resampleSignal(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // Check input and output arguments
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("resampler_mex:invalidInputs", "Three input arguments required.");
    }

    if (nlhs != 1) {
        mexErrMsgIdAndTxt("resampler_mex:invalidOutputs", "One output argument required.");
    }

    // Get the input signal, conversion rate, and interpolation order from MATLAB
    double* inputSignal = mxGetPr(prhs[0]);
    int inputLength = static_cast<int>(mxGetNumberOfElements(prhs[0]));
    double conversionRate = mxGetScalar(prhs[1]);
    int order = static_cast<int>(mxGetScalar(prhs[2]));

    // Create the Resampler object
    Resampler resampler;

    // Calculate the output signal length
    int outputLength = static_cast<int>(ceil(inputLength / conversionRate));

    // Create an output array for the resampled signal
    mxArray* outputSignalArray = mxCreateDoubleMatrix(1, outputLength, mxREAL);
    double* outputSignal = mxGetPr(outputSignalArray);

    // Perform resampling
    resampler.resample(inputSignal, outputSignal, inputLength, conversionRate, order);

    // Return the resampled signal as the output
    plhs[0] = outputSignalArray;
}

// mexFunction entry point
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // Check for the correct number of input and output arguments
    if (nlhs > 1 || nrhs < 3 || nrhs > 3) {
        mexErrMsgIdAndTxt("resampler_mex:usage", "Usage: outputSignal = resampler_mex(inputSignal, conversionRate, order)");
    }

    // Call the resampleSignal function to perform resampling
    resampleSignal(nlhs, plhs, nrhs, prhs);
}
