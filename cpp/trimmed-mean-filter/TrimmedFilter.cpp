// A filter for moving averaging, median, trimmed mean, and weighted median filtering

/*
     Open Source ECG Toolbox, version 2.0, July 2008
 	 Released under the GNU General Public License
 	 Copyright (C) 2008  Reza Sameni
 	 Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, INPG, Grenoble, France
 	 reza.sameni@gmail.com
 	 
 	 This program is free software; you can redistribute it and/or modify it
 	 under the terms of the GNU General Public License as published by the
 	 Free Software Foundation; either version 2 of the License, or (at your
 	 option) any later version.
 	 This program is distributed in the hope that it will be useful, but
 	 WITHOUT ANY WARRANTY; without even the implied warranty of
 	 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 	 Public License for more details.
*/

#include<stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
using namespace std;

#include "list.h"

///////////////////////////////////////////////////////////////////////////////////////////
int TrimmedFilter(const double *arr, double *output, int n, char* type, int w, int alpha=0, const double *h=0)
{
	int i;
	
    // Moving average filter
	if(!strcmp(type,"mean"))
	{
		double temp;
		int n1,n2,i,j;
		for (i=0;i<w/2;i++)
		{
			temp = 0;
            n1 = 0;
			n2 = (i < n-w/2) ? i+(w+1)/2:n;
			for(j=n1;j<n2;j++)
				temp += arr[j];
			output[i] = temp/(n2-n1);
		}
		for (i=w/2;i<n-w/2;i++)
		{
			temp = 0;
			n1 = i-w/2;
			n2 = i+(w+1)/2;
			for(j=n1;j<n2;j++)
				temp += arr[j];
			output[i] = temp/(n2-n1);
		}
		for (i=n-w/2;i<n;i++)
		{
			temp = 0;
			n1 = i-w/2;
			n2 = n;
			for(j=n1;j<n2;j++)
				temp += arr[j];
			output[i] = temp/(n2-n1);
		}
	}
    // Median filter
	else if(!strcmp(type,"median"))
	{
		list lst(w);
		lst.initialize(arr[0]);
		if(w%2==0)
		{
			for(i=0;i<n-(w-1)/2;i++)
			{
				lst.insert(arr[i+(w-1)/2]);
				output[i] = lst.medianEven();
			}
			for(i=n-(w-1)/2;i<n;i++)
			{
				lst.insert(arr[n-1]);
				output[i] = lst.medianEven();
			}
		}
		else
		{
			for(i=0;i<n-(w-1)/2;i++)
			{
				lst.insert(arr[i+(w-1)/2]);
				output[i] = lst.medianOdd();
			}
			for(i=n-(w-1)/2;i<n;i++)
			{
				lst.insert(arr[n-1]);
				output[i] = lst.medianOdd();
			}
		}
	}
    // Trimmed-mean filter
	else if(!strcmp(type,"trmean"))
	{
		list lst(w);
		lst.initialize(arr[0]);
		for(i=0;i<n-(w-1)/2;i++)
		{
			lst.insert(arr[i+(w-1)/2]);
			output[i] = lst.mean(alpha);
		}
		for(i=n-(w-1)/2;i<n;i++)
		{
			lst.insert(arr[n-1]);
			output[i] = lst.mean(alpha);
		}
	}
    // Weighted median filter
	else if(!strcmp(type,"wmedian"))
	{
		list lst(w);
		lst.initialize(arr[0]);
        for(i=0;i<n-(w-1)/2;i++)
        {
            lst.insert(arr[i+(w-1)/2]);
            output[i] = lst.wmedian(h,w,alpha);
        }
        for(i=n-(w-1)/2;i<n;i++)
        {
            lst.insert(arr[n-1]);
            output[i] = lst.wmedian(h,w,alpha);
        }
    }
    else
    {
        return 1;
    }
    
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
// For test
/*
void main(void)
{
double x[14] = {5.7,7.5, 8, 8, 9, 17. , -1 , 22 , 32 , 1.5 , 2 , 9, -18, 7};
//double x[14] = {5.7, 7.5, 8, 8.5, 9, 17, 22, 32, 33, 40, 49, 55, 60, 65};
int n = 14;
double y[14];
int i;

  TrimmedFilter(x,y,n,"trmean",5,2);
  //TrimmedFilter(x,y,n,"mean",3);
  //TrimmedFilter(x,y,n,"median",3);
  
	for(i=0;i<n;i++)
	{
	cout<<i<<'\t'<<x[i]<<'\t'<<y[i]<<endl;
	
	  }
	  }
*/
///////////////////////////////////////////////////////////////////////////////////////////
#include "mex.h"
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	int n,w,a;
	double *x,*y,*h;
	char* type;
	
	if(nrhs < 4)
		mexErrMsgTxt("At least 4 input parameters required.");
	else if(nrhs > 6)
		mexErrMsgTxt("Too many input parameters.");
	
	x = (double *)mxGetPr(prhs[0]);
	n = int(*mxGetPr(prhs[1]));
	type = mxArrayToString(prhs[2]);

	w = int(*mxGetPr(prhs[3]));
	if(w>n)
		mexErrMsgTxt("Inappropriate window length.");

    if(nrhs>4)
		a = int(*mxGetPr(prhs[4]));
	else
		a = 0;

    if(nrhs>5)
		h = (double *)mxGetPr(prhs[5]);
	else
		h = 0;

    if(!strcmp(type,"wmedian") && h==0)
        mexErrMsgTxt("Weighting window required.");
	
	plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
	y = (double *)mxGetPr(plhs[0]);
	
    if(TrimmedFilter(x,y,n,type,w,a,h))
        mexErrMsgTxt("Unknown filter type.");
        
}

///////////////////////////////////////////////////////////////////////////////////////////

