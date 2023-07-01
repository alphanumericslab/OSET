// A class of nonlinear filters using median, trimmed mean, and weighted median filters

#include<stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
using namespace std;

#include "list.h"
#include "nonlinear_filter.h"

///////////////////////////////////////////////////////////////////////////////////////////
NonlinearFilter::NonlinearFilter(const char* ttype, const int ww, const int aalpha, const int bbeta, const double *hh)
{

	int i;
	type = new char[strlen(ttype)+1]; // length including '\0'
	strcpy(type,ttype);

	w = ww;
	alpha = aalpha;
	beta = bbeta;

	lst = new list(w);


	if(!strcmp(type,"wmedian"))
	{
		h = new double[w];
		for(i=0; i<w; i++)
		{
			h[i] = hh[i];
		}
	}
	else
	{
		h = 0;
	}

	this->reset();
}

///////////////////////////////////////////////////////////////////////////////////////////
NonlinearFilter::~NonlinearFilter()
{
	delete type;
	delete h;
	delete lst;
}

///////////////////////////////////////////////////////////////////////////////////////////
void NonlinearFilter::reset(void)
{
	initializedflag = 0; // initial values no longer valid
}
///////////////////////////////////////////////////////////////////////////////////////////
int NonlinearFilter::process(const double *arr, double *output, const int n)
{
	int i;
	if(initializedflag==0)// list not initialized yet
	{
		lst->initialize(arr[0]);
		initializedflag = 1; // list initialized
	}
	
	if(!strcmp(type,"median")) // Median filter
	{
		if(w%2==0)
		{
			for(i=0;i<n;i++)
			{
				lst->insert(arr[i]);
				output[i] = lst->medianEven();
			}
		}
		else
		{
			for(i=0;i<n;i++)
			{
				lst->insert(arr[i]);
				output[i] = lst->medianOdd();
			}
		}
	}
	else if(!strcmp(type,"trmean")) // Trimmed-mean filter
	{
		for(i=0;i<n;i++)
		{
			lst->insert(arr[i]);
			output[i] = lst->mean(alpha,beta);
		}
	}
	else if(!strcmp(type,"wmedian")) // Weighted median filter
	{
        for(i=0;i<n;i++)
        {
            lst->insert(arr[i]);
            output[i] = lst->wmedian(h,w,alpha,beta);
        }
    }
    else
    {
        return 1; // filter type not valid
    }
    
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
/*
// Test bench (from file or command line)
int main(int argc, char* argv[])
{
	int i, cntr, L, w, alpha, beta;
	char *type;
	double d;
	double *x, *y;
	double *h;
	
	//	ifstream fin(argv[1]);
	ifstream fin("testinput.txt");
	if(!fin)
	{
		cerr<<"can not open input file...."<<endl;
		return 1;
	}
	
	//	ifstream fout(argv[2]);
	ofstream fout("testoutput.txt");
	fout.precision(18);
	if(!fout)
	{
		cerr<<"can not open output file...."<<endl;
		return 1;
	}
	
	// type = new char[strlen((char *)argv[3])];
	// strcpy(type,(char *)argv[3]);

	char str[] = "trmean"; // "median", "wmedian", 
	type = new char[strlen(str)];
	strcpy(type,str);

	// L = (int)argv[4];
	L = 750;

	// int w = (int)argv[5];
	w = 50;

	if(w>L)
	{
		cerr<<"data length should be more than filter length...."<<endl;
		return 1;
	}

	//if(!strcmp(type,"trmean") && argc<8)
	//{
	//	cerr<<"trimmed length not specified...."<<endl;
	//	return 1;
	//}

	//alpha = (int)argv[6];
	alpha = 10;
	
	//beta = (int)argv[7];
	beta = 10;

	//	ifstream fcoefs(argv[8]);
	ifstream fcoefs("filtercoefs.txt");
	if(!fcoefs)
	{
		cerr<<"can not open filter coefficient file...."<<endl;
		return 1;
	}

	h = new double[w];
	for(i=0; i<w; i++)
	{
		fcoefs>>h[i];
	}
	

	if( (alpha+beta) > w)
	{
		cerr<<"trimmed length should be less than filter length...."<<endl;
		return 1;
	}

	NonlinearFilter filter(type, w, alpha, beta, h);
	
	
	
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
#include "mex.h"
#define MAX(a,b) ((a)>(b) ? (a):(b))
#define MIN(a,b) ((a)<(b) ? (a):(b))
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	int n,w,a,b;
	double *x,*y,*h;
	char* type;
	NonlinearFilter* filter;
	
	if(nrhs < 5)
	{
		mexErrMsgTxt("At least 5 input parameters required.");
	}
	else if(nrhs > 6)
	{
		mexErrMsgTxt("Too many input parameters.\n");
	}
	
	if (mxGetN(prhs[0]) != 1 && mxGetM(prhs[0]) !=1)
	{
		mexErrMsgTxt("x must be a vector.\n");
	}
	
	if (mxGetN(prhs[2]) != 1 || mxGetM(prhs[2]) !=1)
	{
		mexErrMsgTxt("w must be a scaler.\n");
	}
	
	if (mxGetN(prhs[3]) != 1 || mxGetM(prhs[3]) !=1)
	{
		mexErrMsgTxt("alpha must be a scaler.\n");
	}

	if (mxGetN(prhs[4]) != 1 || mxGetM(prhs[4]) !=1)
	{
		mexErrMsgTxt("beta must be a scaler.\n");
	}

	x = (double *)mxGetPr(prhs[0]);
	n = MAX(mxGetM(prhs[0]), mxGetN(prhs[0]));
	type = mxArrayToString(prhs[1]);
	w = mxGetScalar(prhs[2]);
	a = mxGetScalar(prhs[3]);
	b = mxGetScalar(prhs[4]);

	if(w > n)
	{
		mexErrMsgTxt("Window length longer than data length.\n");
	}
	
	if(!strcmp(type,"wmedian") && nrhs<6)
	{
		mexErrMsgTxt("Weighting window required for 'wmedian' filter type.\n");
	}

	if(!strcmp(type,"wmedian"))
	{
		if( MAX(mxGetM(prhs[5]),mxGetN(prhs[5])) != w)
		{
			mexErrMsgTxt("Inappropriate window length.\n");
		}
		else
		{
			h = (double *)mxGetPr(prhs[5]);
		}
	}
	else
	{
		h = 0;
	}
	
	
	plhs[0] = mxCreateDoubleMatrix(1,n,mxREAL);
	y = (double *)mxGetPr(plhs[0]);

	filter = new NonlinearFilter(type, w, a, b, h);

	if(filter->process(x, y, n))
	{
	  mexErrMsgTxt("Unknown filter type.");
	}
}
///////////////////////////////////////////////////////////////////////////////////////////
