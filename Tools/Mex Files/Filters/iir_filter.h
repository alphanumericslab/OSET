#include <cstdlib>

#include<fstream.h>

#ifndef IIR_FILTER_H_
#define IIR_FILTER_H_

#include "filter.h"

///////////////////////////////////////////////////////////////////////////////////////////
class IIR_filter : public Filter
{
private:
	double *b, *a, *pre_x, *pre_y;
	int M, N;
protected:
	
public:
	//constructor
	//the vectors b,a with length M,N contain the FIR or IIR filter coeficients 
	IIR_filter();
	IIR_filter(const double *b_arg, const int M_arg, const double *a_arg, const int N_arg);
	//destructor
	~IIR_filter();
	//A function for resetting the initial values of the filters
	void reset(void);
	//A function for getting the initial values of the filters
	void get_init(const double *x, const double *y);
	//A function,which filters the data vector x, given its length L
	//y is output vector
	int process(const double *x, double *y, const int L);
};
///////////////////////////////////////////////////////////////////////////////////////////

#endif /*IIR_FILTER_H_*/
