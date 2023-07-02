// Resampler.h: interface for the Resampler class.
// Reza Sameni
// reza.sameni@gmail.com
//
//////////////////////////////////////////////////////////////////////

#if !defined(RESAMPLER_H_)
#define RESAMPLER_H_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include <math.h>
#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>


class Resampler  
{
	int *PreviousSamples;
public:
	
	int volume_coef;
	
	double phase;
				// phase is the fraction t(sec)/Ts(in)
	int Resample(double* x,double* y,int InputLength,double ConversionRate,int Order=2);
	int Resample(int* x,short* y,int InputLength,double ConversionRate,int Order=2);
				// ConversionRate is the fraction Ts(out)/Ts(in) or Fs(in)/Fs(out)
	Resampler();
	virtual ~Resampler();

};

#endif // !defined(RESAMPLER_H_)
