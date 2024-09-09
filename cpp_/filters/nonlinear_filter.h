#include "filter.h"

///////////////////////////////////////////////////////////////////////////////////////////
class NonlinearFilter : public Filter
{
private:
	double *h;
	char* type;
	list *lst;
	int w, alpha, beta;
	int initializedflag;

protected:
	
public:
	//constructor
	//the vectors b,a with length M,N contain the FIR or IIR filter coeficients 
	NonlinearFilter();
	NonlinearFilter(const char* type, const int w, const int alpha=0, const int beta=0, const double *h=0);

	//destructor
	~NonlinearFilter();

	//A function for resetting the initial values of the filters
	void reset(void);
	
	//A function for getting the initial values of the filters
	void get_init(const double *x, const double *y);
	
	//A function, which filters the data vector x, given its length L
	//y is output vector
	int process(const double *x, double *y, const int L);
};
