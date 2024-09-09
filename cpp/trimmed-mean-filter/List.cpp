// An ascending sorted list class used for efficient max/min finding and trimmed mean/median filtering

/*
% Revision History:
%   2008: First release
%   2023: Renamed from deprecated version List.cpp
%   2024: Added max/min finding
%
% Reza Sameni, 2022-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

*/

#include<stdio.h>
//#include <fstream>
//#include <iostream>
//using namespace std;

#include "list.h"

///////////////////////////////////////////////////////////////////////////////////////////
list::list(int N, double init, int init_index)
{
	list_len = N;
	itmlist = new item[N];
	write_pointer = 0;
	head_pointer = 0;
	tail_pointer = list_len-1;
	int i;
	for(i=0;i<list_len;i++)
	{
		itmlist[i].value = init;
		itmlist[i].value_index = init_index;
		itmlist[i].from = i-1;
		itmlist[i].to = i+1;
	}
	itmlist[list_len-1].to = -1;
}

///////////////////////////////////////////////////////////////////////////////////////////
list::~list()
{
	delete[] itmlist;
}

///////////////////////////////////////////////////////////////////////////////////////////
int list::initialize(double init, int init_index)
{
	write_pointer = 0;
	head_pointer = 0;
	tail_pointer = list_len - 1;
	int i;
	for(i=0;i<list_len;i++)
	{
		itmlist[i].value = init;
		itmlist[i].value_index = init_index;
		itmlist[i].from = i-1;
		itmlist[i].to = i+1;
	}
	itmlist[list_len-1].to = -1;
	
	/*	cout<<'\n'<<write_pointer<<'\n'<<"index"<<'\t'<<"value"<<'\t'<<"from"<<'\t'<<"to"<<'\t'<<"head_pointer"<<'\t'<<"tail_pointer"<<endl;
	for(i=0;i<list_len;i++)
	{
	cout<<i<<'\t'<<itmlist[i].value<<'\t'<<itmlist[i].from<<'\t'<<itmlist[i].to<<'\t'<<head_pointer<<'\t'<<tail_pointer<<endl;
	}
	*/	
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
int list::insert(double data)
{
	itmlist[write_pointer].value = data;
	this->resort();
	write_pointer = (write_pointer+1)%list_len;
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
int list::insert(double data, int data_index)
{
	itmlist[write_pointer].value = data;
	itmlist[write_pointer].value_index = data_index;
	this->resort();
	write_pointer = (write_pointer+1)%list_len;
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
int list::resort(void)
{
	int i,temp;
	
	// Moving forward in the list
	while(itmlist[write_pointer].to > -1 && itmlist[write_pointer].value > itmlist[itmlist[write_pointer].to].value)
	{
		temp = itmlist[write_pointer].to;
		if(itmlist[temp].to > -1)
		{
			itmlist[itmlist[temp].to].from = write_pointer;
		}
		itmlist[write_pointer].to = itmlist[temp].to;
		itmlist[temp].from = itmlist[write_pointer].from;
		if(itmlist[write_pointer].from > -1)
		{
			itmlist[itmlist[write_pointer].from].to = temp;
		}
		itmlist[write_pointer].from = temp;
		itmlist[temp].to = write_pointer;
	}
	// Moving backward in the list
	while(itmlist[write_pointer].from > -1 && itmlist[write_pointer].value < itmlist[itmlist[write_pointer].from].value)
	{
		temp = itmlist[write_pointer].from;
		if(itmlist[temp].from > -1)
		{
			itmlist[itmlist[temp].from].to = write_pointer;
		}
		itmlist[write_pointer].from = itmlist[temp].from;
		itmlist[temp].to = itmlist[write_pointer].to;
		if(itmlist[write_pointer].to > -1)
		{
			itmlist[itmlist[write_pointer].to].from = temp;
		}
		itmlist[write_pointer].to = temp;
		itmlist[temp].from = write_pointer;
	}
	
	for(i=0 ; i < list_len ; i++)
	{
		if(itmlist[i].from==-1)
		{
			head_pointer = i;
			break;
		}
	}
	
	for(i=0 ; i< list_len ; i++)
	{
		if(itmlist[i].to==-1)
		{
			tail_pointer = i;
			break;
		}
	}
	
	
	/*	cout<<'\n'<<write_pointer<<'\n'<<"index"<<'\t'<<"value"<<'\t'<<"from"<<'\t'<<"to"<<'\t'<<"head_pointer"<<'\t'<<"tail_pointer"<<endl;
	for(i=0;i<list_len;i++)
	{
	cout<<i<<'\t'<<itmlist[i].value<<'\t'<<itmlist[i].from<<'\t'<<itmlist[i].to<<'\t'<<head_pointer<<'\t'<<tail_pointer<<endl;
	}
	*/	
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
double list::mean(int a) // neglects the first and last 'a' samples 
{
	int position_pointer = head_pointer;
	int cntr = 0;
	double sum = 0;
	while(cntr < list_len-a){
		if(cntr >= a){
			sum += itmlist[position_pointer].value;
		}
		cntr++;
		position_pointer = itmlist[position_pointer].to;
	}

	return sum/(cntr-a);
}

///////////////////////////////////////////////////////////////////////////////////////////
double list::wmedian(const double *weight, int m, int a) // neglects the first and last 'a' samples 
{
	int position_pointer = head_pointer;
	int cntr = 0;
	double sum = 0;
	while(cntr < list_len-a){
		if(cntr >= a){
			sum += weight[cntr]*itmlist[position_pointer].value;
		}
		cntr++;
		position_pointer = itmlist[position_pointer].to;
	}

	return sum;
}

///////////////////////////////////////////////////////////////////////////////////////////
double list::median_odd(void)
{
	int position_pointer = head_pointer;
	int cntr = 0;
	while(cntr < list_len/2)
	{
		position_pointer = itmlist[position_pointer].to;
		cntr++;
	}
	return itmlist[position_pointer].value;
}

///////////////////////////////////////////////////////////////////////////////////////////
double list::median_even(void)
{
	int position_pointer = head_pointer;
	int cntr = 0;
	while(cntr < list_len/2)
	{
		position_pointer = itmlist[position_pointer].to;
		cntr++;
	}
	return (itmlist[position_pointer].value + itmlist[itmlist[position_pointer].from].value)/2.;
}

///////////////////////////////////////////////////////////////////////////////////////////
double list::max(void)
{
	return itmlist[tail_pointer].value;
}

///////////////////////////////////////////////////////////////////////////////////////////
int list::max_index(void)
{
	return itmlist[tail_pointer].value_index;
}

///////////////////////////////////////////////////////////////////////////////////////////
double list::min(void)
{
	return itmlist[head_pointer].value;
}

///////////////////////////////////////////////////////////////////////////////////////////
int list::min_index(void)
{
	return itmlist[head_pointer].value_index;
}
