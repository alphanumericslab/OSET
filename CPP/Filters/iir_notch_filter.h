#ifndef IIR_NOTCH_FILTER_H_
#define IIR_NOTCH_FILTER_H_

//#include <gsl/gsl_math.h>
#define M_PI 4*atan(1)
#include <math.h>

#include <cstdlib>
#include<fstream.h>


#include "iir_filter.h"

/////////////////////////////////////////////////////////////////////////
class IIR_notch_filter: public IIR_filter
{
private:
	double *a_cof, *b_cof;
//	double alpha, w0;
public:
	IIR_notch_filter(): IIR_filter()
	{}
	//A constructor for setting the notch frequency f0 and the cutoff frequency fc
	IIR_notch_filter(const double f0, const double fc, const double fs = 1000);
	~IIR_notch_filter();
};
/////////////////////////////////////////////////////////////////////////

#endif /*IIR_NOTCH_FILTER_H_*/
