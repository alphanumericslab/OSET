// A single or two-stage median filter

/*
     Open Source ECG Toolbox, version 2.0, March 2008
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

double median(double arr [],int n);
void median_filter_batch(const double *arr,double *output,int n,int w);

///////////////////////////////////////////////////////////////////////////////////////////
#include "mex.h"
void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	int n,w1,w2;
	double *x,*y,*ytemp;

	if(nrhs == 0)
	{
		mexPrintf(" y = Median(x,N,L1,L2),\n Single or Double stage median filter\n inputs:\n\t x: input vector. x should be in 'double' format.\n\t N: input signal length (length(x))\n\t L1: first stage median window length (in samples)\n\t L2: second stage median window length (in samples)\n\t approach:\n output:\n\t y: output column vector\n\n\t Open Source ECG Toolbox, version 2.0, March 2008\n\t Released under the GNU General Public License\n\t Copyright (C) 2008  Reza Sameni\n\t Sharif University of Technology, Tehran, Iran -- GIPSA-Lab, INPG, Grenoble, France\n\t reza.sameni@gmail.com\n\t \n\t This program is free software; you can redistribute it and/or modify it\n\t under the terms of the GNU General Public License as published by the\n\t Free Software Foundation; either version 2 of the License, or (at your\n\t option) any later version.\n\t This program is distributed in the hope that it will be useful, but\n\t WITHOUT ANY WARRANTY; without even the implied warranty of\n\t MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General\n\t Public License for more details.\n");
	}
	else{
		if(nrhs < 3)
		{
			mexErrMsgTxt("At least 3 input parameters required.");
		}
		else if(nrhs > 4)
		{
			mexErrMsgTxt("Too many input parameters.");
		}
		
		n = int(*mxGetPr(prhs[1]));
		
		
		x = (double *)mxGetPr(prhs[0]);
		
		
		if(nrhs==3)	// single pass median filter
		{
			w1 = int(*mxGetPr(prhs[2]));
			plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
			y = (double *)mxGetPr(plhs[0]);
			
			median_filter_batch(x,y,n,w1);
		}
		else if(nrhs==4) // double pass median filter
		{
			w1 = int(*mxGetPr(prhs[2]));
			w2 = int(*mxGetPr(prhs[3]));

			plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
			y = (double *)mxGetPr(plhs[0]);

			ytemp = (double *)mxGetPr(mxCreateDoubleMatrix(n,1,mxREAL));
			
			median_filter_batch(x,ytemp,n,w1);
			median_filter_batch(ytemp,y,n,w2);
		}
	}
}
///////////////////////////////////////////////////////////////////////////////////////////
/*
#include<stdio.h>
#include <fstream>
#include <iostream>
using namespace std;
void main(void)
{
double x[14] = {5.7,7.5, 8, 8, 9, 17. , -1 , 22 , 32 , 1.5 , 2 , 9, -18, 7};
int i,n = 14;
double y[14];

  median_filter_batch(x,y,n,5);
  	for(i=0;i<n;i++)
  	{
  		cout<<i<<'\t'<<x[i]<<'\t'<<y[i]<<endl;
  	}
  }
  
	/*int main(int argc,char *argv[])
	{
	if(argc==5)
	{
	median_filter_batch((double*)argv[1],(double*)argv[2],(int)argv[3],(int)argv[4]);
	return 1;
	}
	else
	return 0;
	}
*/
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////
#define SWAP(a,b) temp = (a);(a) = (b);(b) = temp; 
double select(int k, int n, double *arr) 
{ 
	int i,ir,j,l,mid; 
	double a,temp;
	l = 0;
	ir = n-1; 
	for(;;){ 
		if(ir <= l+1){									// Active partition contains 1 or 2 elements, 
			if (ir == l+1 && arr[ir] < arr[l]) {		// Case of 2 elements. 
				SWAP(arr[l],arr[ir])
			}
			return arr[k]; 
		} else{
			mid = (l + ir) >> 1;						// Choose median of left, center, and right el- 
			SWAP(arr[mid],arr[l+1])						// ements as partitioning element a. Also 
				if(arr[l] > arr[ir]){					// rearrange so that arr[l] < arr [1+1], 
					SWAP(arr[l],arr[ir])				// arr[ir] > arr[l+1]. 
				} 
				if(arr [l+1] > arr[ir]){ 
					SWAP (arr[l+1],arr[ir]) 
				} 
				if(arr[l] > arr[l+1]){ 
					SWAP(arr [l],arr[l+1]) 
				} 
				i = l+1;								// Initialize pointers for partitioning. 
				j = ir;
				a = arr[l+1];							// Partitioning element. 
				for(;;){								// Beginning of innermost loop. 
					do i++; while (arr[i]<a);			// Scan up to find element > a. 
					do j--; while (arr[j]>a);			// Scan down to find element < a. 
					if(j<i) break;						// Pointers crossed. Partitioning complete. 
					SWAP(arr[i],arr[j]);
				}										// End of innermost loop. 
				arr[l+1] = arr[j];						// Insert partitioning element. 
				arr[j] = a;
				if(j>= k) ir= j-1;						// Keep active the partition that contains the kth element. 
				if(j<= k) l = i;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
double median(double *arr,int n)
{
	if(n%2==0)
	{
		return (select(n/2-1,n,arr) + select(n/2,n,arr))/2.;
	}
	else
	{
		return select((n-1)/2,n,arr);
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
void median_filter_batch(const double *arr,double *output,int n,int w)
{
	double *temparray;
	temparray = new double[w];
	int start,stop,i,j;
	// cout<<"i"<<'\t'<<"start"<<'\t'<<"stop"<<'\t'<<"arr[i]"<<'\t'<<"out[i]"<<'\t'<<"stop-start"<<endl;
	for (i=0;i<n;i++)
	{
		start = (i > w/2) ? i-w/2:0;
		stop = (i < n-w/2) ? i+(w+1)/2:n;
		for(j=start;j<stop;j++)
			temparray[j-start] = arr[j];
		output[i] = median(temparray,stop-start);
		// cout<<i<<'\t'<<start<<'\t'<<stop<<'\t'<<arr[i]<<'\t'<<output[i]<<'\t'<<stop-start<<endl;
	}
	delete[] temparray;
}
///////////////////////////////////////////////////////////////////////////////////////////
