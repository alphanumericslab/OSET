// Resampler.h: interface for the Resampler class.

/*
% Revision History:
%   2012: First release
%   2023: Renamed from deprecated version Resampler.cpp
%
% Reza Sameni, 2012-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

*/
//////////////////////////////////////////////////////////////////////

#if !defined(RESAMPLER_H_)
#define RESAMPLER_H_

#include <iostream>
#include <fstream>
#include <stdlib.h>


class Resampler  
{
	int *previous_samples;
public:
		
	double phase;
				// phase is the fraction t(sec)/Ts(in)
	int resample(double* x, double* y, int InputLength, double ConversionRate, int Order=2);
	int resample(int* x, int* y, int InputLength, double ConversionRate, int Order=2);
				// ConversionRate is the fraction Ts(out)/Ts(in) or Fs(in)/Fs(out)
	Resampler();
	virtual ~Resampler();

};

#endif // !defined(RESAMPLER_H_)
