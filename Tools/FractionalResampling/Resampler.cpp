// Resampler.cpp: implementation of a fractional Resampler class.
// Reza Sameni
// reza.sameni@gmail.com
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Resampler.h"
#include <math.h>


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


Resampler::~Resampler()
{
	delete[] PreviousSamples;
}

Resampler::Resampler()
{
	PreviousSamples = new int[4];
	for(int i=0;i<4;i++)
		PreviousSamples[i] = 0;
	phase=0;
}
//////////////////////////////////////////////////////////////////////

int Resampler::Resample(double* x,double* y,int InputLength,double ConversionRate,int Order)
{
	int k,i;
	double y1,y2,y3,y4;
	int OutputLength = (int)ceil(InputLength / ConversionRate);;

	switch(Order){

	case 2:// Second Order Resampler (Triangular Convolution)
		{
			for(i=0;i<OutputLength;i++){

				k=(int)ceil(phase-.5);
				y[i]=(k+1-phase)*x[k] + (phase-k)*x[k+1];
				phase += ConversionRate;
			}
			break;
		}
	case 3:// Third Order Resampler (Quadratic Convolution)
		{
			for(i=0;i<OutputLength;i++){
				phase += ConversionRate;
				k=(int)ceil(phase - 1.5);

				y1=x[k];
				y2=x[k+1];
				y3=x[k+2];

				y[i]=(phase-k)/2 * ((y1-2*y2+y3)*(phase-k) - (3*y1-4*y2+y3) ) + y1;

			}
			break;
		}
	case 4:// Fourth Order Resampler
		{
			for(i=0;i<OutputLength;i++){
				phase += ConversionRate;
				k=(int)ceil(phase - 2.);

				y1=x[k];
				y2=x[k+1];
				y3=x[k+2];
				y4=x[k+3];

				y[i]=- (y1-3*y2+3*y3-y4)*(phase-k)*(phase-k)*(phase-k)/6
					+ (2*y1-5*y2+4*y3-y4)*(phase-k)*(phase-k)/2
					- (11*y1-18*y2+9*y3-2*y4)*(phase-k)/6
					+ y1;

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
int Resampler::Resample(int* x,short* y,int InputLength,double ConversionRate,int Order)
{
	int k,i;
	int y1,y2,y3,y4;
	int OutputLength = (int)floor(InputLength / ConversionRate) -1;

	switch(Order){

	case 2:// Second Order Resampler (Triangular Convolution)
		{
			phase -= floor(phase);
			for(i=0;i<OutputLength;i++){

				k=(int)ceil(phase-.5);
				if(k<InputLength)
					y[i]=(short)floor(.5 + (k+1-phase)*x[k] + (phase-k)*x[k+1]);
				else
					y[i]=(short)floor(.5 + (k+1-phase)*x[k] + (phase-k)*0);

				phase += ConversionRate;
			}
			break;
		}
	case 3:// Third Order Resampler (Quadratic Convolution)
		{
			phase-=floor(phase);
			for(i=0;i<OutputLength;i++){
				phase += ConversionRate;
				k=(int)ceil(phase - 1.5);

				y1=(k<InputLength && k>0) ? x[k] : 0;
				y2=((k+1)<InputLength && k>0) ? x[k+1] : 0;
				y3=((k+2)<InputLength && k>0) ? x[k+2] : 0;

				y[i]=(short)floor(.5 + (phase-k)/2. * ((y1-2*y2+y3)*(phase-k) - (3*y1-4*y2+y3) ) + y1);

			}
			break;
		}
	case 4:// Fourth Order Resampler
		{
			phase -= floor(phase);
			for(i=0;i<OutputLength;i++){
				phase += ConversionRate;
				k=(int)ceil(phase - 2.);

				y1 = (k<InputLength) ? ((k>-1) ? x[k] : PreviousSamples[3 + k]) : 0;
				y2 = ((k+1)<InputLength) ? ((k>-2) ? x[k+1] : PreviousSamples[4 + k]) : 0;
				y3 = ((k+2)<InputLength) ? ((k>-3) ? x[k+2] : PreviousSamples[5 + k]) : 0;
				y4 = ((k+3)<InputLength) ? ((k>-4) ? x[k+3] : PreviousSamples[6 + k]) : 0;

				y[i]=(short)floor( (1/float(volume_coef))*(.5 - (y1-3*y2+3*y3-y4)*(phase-k)*(phase-k)*(phase-k)/6.
					+ (2*y1-5*y2+4*y3-y4)*(phase-k)*(phase-k)/2.
					- (11*y1-18*y2+9*y3-2*y4)*(phase-k)/6.
					+ y1));

			}
			break;
		}
	default:
		{
			break;
		}
	}
	
	for(i=0;i<4;i++)
		PreviousSamples[i] = x[InputLength-4+i];

	return 1;
}
