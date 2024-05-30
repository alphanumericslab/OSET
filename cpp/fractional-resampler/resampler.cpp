// Resampler.cpp: implementation of a fractional Resampler class.

/*
% Revision History:
%   2012: First release
%   2023: Renamed from deprecated version Resampler.cpp
%
% Reza Sameni, 2012-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

*/
//////////////////////////////////////////////////////////////////////

#include "resampler.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

// Destructor
Resampler::~Resampler()
{
    // Clean up dynamically allocated memory
    delete[] previous_samples;
}

// Constructor
Resampler::Resampler()
{
    // Allocate memory for previous_samples array to store past samples
    previous_samples = new int[4];

    // Initialize previous_samples array with zeros
    for (int i = 0; i < 4; i++)
        previous_samples[i] = 0;

    // Initialize the phase to zero
    phase = 0;
}

//////////////////////////////////////////////////////////////////////

// Function for resampling a continuous double precision input signal
// x: Input signal
// y: Output signal (resampled)
// InputLength: Length of the input signal 'x'
// ConversionRate: Conversion rate for resampling (input sample rate / output sample rate)
// Order: Order of interpolation (2: Second Order, 3: Third Order, 4: Fourth Order)
int Resampler::resample(double* x, double* y, int InputLength, double ConversionRate, int Order)
{
    int k, i;
    double y1, y2, y3, y4;
    int OutputLength = (int)ceil(InputLength / ConversionRate);

    switch (Order)
    {
        case 2: // Second Order Resampler (Triangular Convolution)
            {
                for (i = 0; i < OutputLength; i++)
                {
                    k = (int)floor(phase);
                    if( (k + 1) < InputLength )
                    {
                        y[i] = (k + 1 - phase) * x[k] + (phase - k) * x[k + 1];
                    }
                    phase += ConversionRate;
                }
                break;
            }
        case 3: // Third Order Resampler (Quadratic Convolution)
            {
                for (i = 0; i < OutputLength; i++)
                {
                    phase += ConversionRate;
                    k = (int)ceil(phase - 1.5);

                    if( (k + 2) < InputLength )
                    {
                        y1 = x[k];
                        y2 = x[k + 1];
                        y3 = x[k + 2];

                        y[i] = (phase - k) / 2 * ((y1 - 2 * y2 + y3) * (phase - k) - (3 * y1 - 4 * y2 + y3)) + y1;
                    }
                }
                break;
            }
        case 4: // Fourth Order Resampler
            {
                for (i = 0; i < OutputLength; i++)
                {
                    phase += ConversionRate;
                    k = (int)ceil(phase - 2.);

                    if( (k + 3) < InputLength )
                    {
                        y1 = x[k];
                        y2 = x[k + 1];
                        y3 = x[k + 2];
                        y4 = x[k + 3];

                        y[i] = -((y1 - 3 * y2 + 3 * y3 - y4) * (phase - k) * (phase - k) * (phase - k)) / 6
                            + (2 * y1 - 5 * y2 + 4 * y3 - y4) * (phase - k) * (phase - k) / 2
                            - (11 * y1 - 18 * y2 + 9 * y3 - 2 * y4) * (phase - k) / 6
                            + y1;
                    }
                }
                break;
            }
        default:
            {
                break;
            }
    }

    return 1;
}

///////////////////////////////////////////////////////////////////////////////
// Function for resampling a continuous integer input signal
// x: Input signal (integer)
// y: Output signal (resampled as integers)
// InputLength: Length of the input signal 'x'
// ConversionRate: Conversion rate for resampling (input sample rate / output sample rate)
// Order: Order of interpolation (2: Second Order, 3: Third Order, 4: Fourth Order)
int Resampler::resample(int* x, int* y, int InputLength, double ConversionRate, int Order)
{
    int k, i;
    int y1, y2, y3, y4;
    int OutputLength = (int)floor(InputLength / ConversionRate) - 1;

    // Adjust phase to ensure proper alignment with input samples
    phase -= floor(phase);
    switch (Order)
    {
        case 2: // Second Order Resampler (Triangular Convolution)
            {

                for (i = 0; i < OutputLength; i++)
                {
                    k = (int)ceil(phase - 0.5);
                    if (k < InputLength)
                        y[i] = (int)floor(0.5 + (k + 1 - phase) * x[k] + (phase - k) * x[k + 1]);
                    else
                        y[i] = (int)floor(0.5 + (k + 1 - phase) * x[k] + (phase - k) * 0);

                    phase += ConversionRate;
                }
                break;
            }
        case 3: // Third Order Resampler (Quadratic Convolution)
            {
                for (i = 0; i < OutputLength; i++)
                {
                    phase += ConversionRate;
                    k = (int)ceil(phase - 1.5);

                    y1 = (k < InputLength && k > 0) ? x[k] : 0;
                    y2 = ((k + 1) < InputLength && k > 0) ? x[k + 1] : 0;
                    y3 = ((k + 2) < InputLength && k > 0) ? x[k + 2] : 0;

                    y[i] = (int)floor(0.5 + (phase - k) / 2. * ((y1 - 2 * y2 + y3) * (phase - k) - (3 * y1 - 4 * y2 + y3)) + y1);
                }
                break;
            }
        case 4: // Fourth Order Resampler
            {
                for (i = 0; i < OutputLength; i++)
                {
                    phase += ConversionRate;
                    k = (int)ceil(phase - 2.);

                    // Handle cases when k is out of bounds and use previous_samples for interpolation
                    y1 = (k < InputLength) ? ((k > -1) ? x[k] : previous_samples[3 + k]) : 0;
                    y2 = ((k + 1) < InputLength) ? ((k > -2) ? x[k + 1] : previous_samples[4 + k]) : 0;
                    y3 = ((k + 2) < InputLength) ? ((k > -3) ? x[k + 2] : previous_samples[5 + k]) : 0;
                    y4 = ((k + 3) < InputLength) ? ((k > -4) ? x[k + 3] : previous_samples[6 + k]) : 0;

                    // Perform interpolation using the selected order of interpolation
                    y[i] = (int)floor(0.5 - (y1 - 3 * y2 + 3 * y3 - y4) * (phase - k) * (phase - k) * (phase - k) / 6.
                                      + (2 * y1 - 5 * y2 + 4 * y3 - y4) * (phase - k) * (phase - k) / 2.
                                      - (11 * y1 - 18 * y2 + 9 * y3 - 2 * y4) * (phase - k) / 6.
                                      + y1);
                }
                break;
            }
        default:
            {
                break;
            }
    }

    // Store the last four samples of the input signal in the previous_samples array for interpolation in the next block
    for (i = 0; i < 4; i++)
        previous_samples[i] = x[InputLength - 4 + i];

    return 1;
}
