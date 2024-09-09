// A filter for moving averaging, median, trimmed mean, and weighted median filtering

/*
% Revision History:
%   2008: First release
%   2023: Renamed from deprecated version TrimmedFilter.cpp
%
% Reza Sameni, 2022-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

*/

#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <cstdlib>

using namespace std;

#include "list.h"
#include "trimmed_filter.h"

/////////////////////
int trimmed_filter(const double *arr, double *output, int n, const char* type, int w, int alpha, const double* h)
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
			for(i = 0 ; i < (w-1)/2 ; i++)
			{
				lst.insert(arr[i]);
			}
			for(i = 0 ; i < n-(w-1)/2 ; i++)
			{
				lst.insert(arr[i+(w-1)/2]);
				output[i] = lst.median_even();
			}
			for(i = n-(w-1)/2 ; i < n ; i++)
			{
				lst.insert(arr[n-1]);
				output[i] = lst.median_even();
			}
		}
		else
		{
			for(i=0;i<(w-1)/2;i++)
			{
				lst.insert(arr[i]);
			}
			for(i=0;i<n-(w-1)/2;i++)
			{
				lst.insert(arr[i+(w-1)/2]);
				output[i] = lst.median_odd();
			}
			for(i=n-(w-1)/2;i<n;i++)
			{
				lst.insert(arr[n-1]);
				output[i] = lst.median_odd();
			}
		}
	}
    // Trimmed-mean filter
	else if(!strcmp(type,"trmean"))
	{
		list lst(w);
		lst.initialize(arr[0]);
		for(i=0;i < (w-1)/2;i++)
		{
			lst.insert(arr[i]);
		}
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
        for(i=0;i < (w-1)/2;i++)
        {
            lst.insert(arr[i]);
        }
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
    // local max finder
	else if(!strcmp(type,"max"))
	{
		list lst(w);
		lst.initialize(arr[0]);
        for(i=0 ; i < (w-1)/2 ; i++)
        {
            lst.insert(arr[i]);
        }
        for(i = 0 ; i < n-(w-1)/2 ; i++)
        {
            lst.insert(arr[i+(w-1)/2]);
            output[i] = lst.max();
        }
        for(i = n-(w-1)/2 ; i < n ; i++)
        {
            lst.insert(arr[n-1]);
            output[i] = lst.max();
        }
    }
    // local min finder
	else if(!strcmp(type,"min"))
	{
		list lst(w);
		lst.initialize(arr[0]);
        for(i = 0 ; i < (w-1)/2 ; i++)
        {
            lst.insert(arr[i]);
        }
        for(i = 0 ; i < n-(w-1)/2 ; i++)
        {
            lst.insert(arr[i+(w-1)/2]);
            output[i] = lst.min();
        }
        for(i = n-(w-1)/2 ; i < n ; i++)
        {
            lst.insert(arr[n-1]);
            output[i] = lst.min();
        }
    }
    else
    {
        return 1;
    }
    
	return 0;
}

/////////////////////
int trimmed_filter_max_min_indexes(const double *arr, int *output, int n, const char* type, int w)
{
	int i;
	
    // local max finder
	if(!strcmp(type,"max_index"))
	{
		list lst(w);
		lst.initialize(arr[0]);
        for(i=0 ; i < (w-1)/2 ; i++)
        {
            lst.insert(arr[i], i);
        }
        for(i = 0 ; i < n-(w-1)/2 ; i++)
        {
            lst.insert(arr[i+(w-1)/2], i+(w-1)/2);
            output[i] = lst.max_index();
        }
        for(i = n-(w-1)/2 ; i < n ; i++)
        {
            lst.insert(arr[n-1], n-1);
            output[i] = lst.max_index();
        }
    }
    // local min finder
	else if(!strcmp(type,"min_index"))
	{
		list lst(w);
		lst.initialize(arr[0]);
        for(i = 0 ; i < (w-1)/2 ; i++)
        {
            lst.insert(arr[i], i);
        }
        for(i = 0 ; i < n-(w-1)/2 ; i++)
        {
            lst.insert(arr[i+(w-1)/2], i+(w-1)/2);
            output[i] = lst.min_index();
        }
        for(i = n-(w-1)/2 ; i < n ; i++)
        {
            lst.insert(arr[n-1], n-1);
            output[i] = lst.min_index();
        }
    }
    else
    {
        return 1;
    }
    
	return 0;
}
