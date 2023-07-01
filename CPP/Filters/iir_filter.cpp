#include "iir_filter.h"

//----------------------------------------------------------------------------------------
IIR_filter::IIR_filter(): b(NULL), a(NULL), M(0), N(0)
{
}
//----------------------------------------------------------------------------------------
IIR_filter::IIR_filter(const double *b_arg, const int M_arg, const double *a_arg, const int N_arg)
{
	b = (double*)b_arg;
	a = (double*)a_arg;
	M = M_arg;
	N = N_arg; 
	pre_x = new double[M];//previous input
	pre_y = new double[N];//previous output
	this->reset();
	//fix me:finding a better way to allocate memory with setting to zero: zare 87/9/3
	//	pre_x = (double*) calloc (M,sizeof(double));//previous input
	//	pre_y = (double*) calloc (N,sizeof(double));//previous output
	//	double pre_x[M],pre_y[N];
}
//----------------------------------------------------------------------------------------
IIR_filter::~IIR_filter()
{
	delete[] pre_x;
	delete[] pre_y;
	//	free (pre_x);
	//	free (pre_y);
}
//----------------------------------------------------------------------------------------
void IIR_filter::reset(void)
{
    for (int i = 0; i<M; i++)
		pre_x[i] = 0;
    for (int j = 0; j<N; j++)
		pre_y[j] = 0;   
}
//----------------------------------------------------------------------------------------
void IIR_filter::get_init(const double *x, const double *y)
{
	x = pre_x;
	y = pre_y;
}
//----------------------------------------------------------------------------------------
int IIR_filter::process(const double *x, double *y, const int L)
{ 
	double inp, out;
	int i, j, n;
	int MN = (M>N) ? M:N;
	
	for (n = 0; n<MN; n++)
    {
        y[n] = 0.0;
        for (i = 0; i<M; i++)  
        {
			inp = (n-i<0) ? pre_x[n-i+M] : x[n-i];	// previous inputs
			y[n] += b[i] * inp;
        }  
        for (j = 1; j<N; j++) 
        {
			out = (n-j<0) ? pre_y[n-j+N] : y[n-j];	// previous outputs
			y[n] -= a[j] * out;
        }
		y[n] /= a[0];
	}
	
	for (n = MN; n<L; n++)
	{
		y[n] = 0.0;
		for (i = 0; i<M; i++)
			y[n] += b[i] * x[n-i];
		
		for (j = 1; j<N; j++)
			y[n] -= a[j] * y[n-j];
		y[n] /= a[0];
	}
	
	// we need the last M elements of x and the last elements of the y vector
	for (i = 0; i<M; i++)
	{
		pre_x[i] = x[L-M+i];
		//cout<<pre_x[i]<<endl;
	}
	//cout<<endl;
	
	for (j = 0; j<N; j++)
		pre_y[j] = y[L-N+j];
	
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
/*
// Test bench (from file or command line)
int main(int argc, char* argv[])
{
	int i, cntr, L, M, N;
	double d;
	double *x, *y;
	double *b, *a;
	
	//	ifstream fcoefs(argv[1]);
	ifstream fcoefs("filtercoefs.txt");
	if(!fcoefs)
	{
		cerr<<"can not open filter coefficient file...."<<endl;
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
	
	fcoefs>>M;
	b = new double[M];
	for(i=0; i<M; i++)
	{
		fcoefs>>b[i];
	}
	
	fcoefs>>N;
	a = new double[N];
	for(i=0; i<N; i++)
	{
		fcoefs>>a[i];
	}
	
	IIR_filter filter(b, M, a, N);
	
	//L = (int)argv[4];
	L = 750;
	if(M>L || N>L)
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
		//cout<<cntr<<endl;
		filter.process(x, y, cntr);
		for(i=0; i<cntr; i++)
		{
			fout<<y[i]<<endl;
		}
	}
	
	fin.close();
	fout.close();
	fcoefs.close();
	
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
	int n, M, N;
	double *x, *y, *b, *a;
	
	if(nrhs < 3 || nrhs > 3)
	{
		mexErrMsgTxt("3 input parameters required.");
	}
	
	if (mxGetN(prhs[0]) != 1 && mxGetM(prhs[0]) !=1)
	{
		mexErrMsgTxt("B must be a vector.\n");
	}
	
	if (mxGetN(prhs[1]) != 1 && mxGetM(prhs[1]) !=1)
	{
		mexErrMsgTxt("A must be a vector.\n");
	}
	
	if (mxGetN(prhs[2]) != 1 && mxGetM(prhs[2]) !=1)
	{
		mexErrMsgTxt("x must be a vector.\n");
	}
	
	b = (double *)mxGetPr(prhs[0]);
	a = (double *)mxGetPr(prhs[1]);
	x = (double *)mxGetPr(prhs[2]);
	
	M = MAX(mxGetM(prhs[0]),mxGetN(prhs[0]));
	N = MAX(mxGetM(prhs[1]),mxGetN(prhs[1]));
	n = MAX(mxGetM(prhs[2]),mxGetN(prhs[2]));
	
	plhs[0] = mxCreateDoubleMatrix(1,n,mxREAL);
	y = (double *)mxGetPr(plhs[0]);
	
	IIR_filter filter(b, M, a, N);
	filter.process(x, y, n);       
}
///////////////////////////////////////////////////////////////////////////////////////////
