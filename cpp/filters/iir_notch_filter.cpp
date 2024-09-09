#include "iir_notch_filter.h"
//----------------------------------------------------------------------------------------
IIR_notch_filter::IIR_notch_filter(const double f0, const double fc, const double fs)
{
	double alpha, w0;
	a_cof = new double[3];
	b_cof = new double[3];
	
	//calculate b,a coefficients
	w0 = 2*M_PI*f0/fs;
	alpha = 1/(1+tan(M_PI*fc/fs));
	
	a_cof[0] = 1;
	a_cof[1] = -2*alpha*cos(w0);
	a_cof[2] = 2*alpha-1;

	b_cof[0] = alpha;
	b_cof[1] = -2*cos(w0)*alpha;
	b_cof[2] = alpha;

	this->IIR_filter::IIR_filter(b_cof, 3, a_cof, 3);//call base constructor
}
//----------------------------------------------------------------------------------------
IIR_notch_filter::~IIR_notch_filter()
{
	delete[] a_cof;
	delete[] b_cof;
}


///////////////////////////////////////////////////////////////////////////////////////////
/*
// Test bench (from file or command line)
int main(int argc, char* argv[])
{
	int i, cntr, L;
	double d, f0, fc, fs;
	double *x, *y;
	
	//	ifstream fcoefs(argv[1]);
	ifstream fparams("notchfilterparams.txt");
	if(!fparams)
	{
		cerr<<"can not open notch filter parameter file...."<<endl;
		return 1;
	}
	
	//	ifstream fin(argv[2]);
	ifstream fin("testinput.txt");
	if(!fin)
	{
		cerr<<"can not open input file...."<<endl;
		return 1;
	}
	
	//	ifstream fout(argv[3]);
	ofstream fout("testoutput.txt");
	fout.precision(18);
	if(!fout)
	{
		cerr<<"can not open output file...."<<endl;
		return 1;
	}
	
	fparams>>f0;
	fparams>>fc;
	fparams>>fs;

	//cout<<f0<<'\t'<<fc<<'\t'<<fs<<endl;

	IIR_notch_filter filter(f0, fc, fs);
	
	//L = (int)argv[4];
	L = 100;

	if(L<3)
	{
		cerr<<"data length should be more than filter length...."<<endl;
		return 1;
	}
	
	
	x = new double[L];
	y = new double[L];
	
	cntr = 0;
	while (!fin.eof())
	{
		if(cntr<L)
		{
			fin>>d;
			x[cntr] = d;
			cntr++;
		}
		else
		{
			filter.process(x, y, L);
			for(i=0; i<L; i++)
			{
				fout<<y[i]<<endl;
			}
			cntr = 0;
		}
	}
	cntr--;
	
	if(cntr>0)
	{
		filter.process(x, y, cntr);
		for(i=0; i<cntr; i++)
		{
			fout<<y[i]<<endl;
		}
	}
	
	fin.close();
	fout.close();
	fparams.close();
	
	return 0;
}
*/
///////////////////////////////////////////////////////////////////////////////////////////
// Matlab mex file
#include "mex.h"
#define MAX(a,b) ((a)>(b) ? (a):(b))
#define MIN(a,b) ((a)<(b) ? (a):(b))
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int n, f0, fc, fs;
	double *x, *y;
	
	if(nrhs < 4 || nrhs > 4)
	{
		mexErrMsgTxt("4 input parameters required.");
	}
	
	if (mxGetN(prhs[0]) != 1 && mxGetM(prhs[0]) !=1)
	{
		mexErrMsgTxt("x must be a vector.\n");
	}
	
	if (mxGetN(prhs[1]) != 1 || mxGetM(prhs[1]) !=1)
	{
		mexErrMsgTxt("f0 must be a scaler.\n");
	}
	
	if (mxGetN(prhs[2]) != 1 || mxGetM(prhs[2]) !=1)
	{
		mexErrMsgTxt("fc must be a scaler.\n");
	}

	if (mxGetN(prhs[3]) != 1 || mxGetM(prhs[3]) !=1)
	{
		mexErrMsgTxt("fs must be a scaler.\n");
	}
	
	x = (double *)mxGetPr(prhs[0]);
	f0 = mxGetScalar(prhs[1]);
	fc = mxGetScalar(prhs[2]);
	fs = mxGetScalar(prhs[3]);
	
	n = MAX(mxGetM(prhs[0]),mxGetN(prhs[0]));
	
	plhs[0] = mxCreateDoubleMatrix(1,n,mxREAL);
	y = (double *)mxGetPr(plhs[0]);
	
	IIR_notch_filter filter(f0, fc, fs);
	filter.process(x, y, n);       
}
///////////////////////////////////////////////////////////////////////////////////////////
