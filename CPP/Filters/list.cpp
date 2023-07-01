#include<stdio.h>
#include <fstream>
#include <iostream>
using namespace std;

#include "list.h"

///////////////////////////////////////////////////////////////////////////////////////////
list::list(int N,double init)
{
	len = N;
	itmlist = new item[N];
	ptr = 0;
	start = 0;
	stop = len-1;
	int i;
	for(i=0;i<len;i++)
	{
		itmlist[i].value = init;
		itmlist[i].from = i-1;
		itmlist[i].to = i+1;
	}
	itmlist[len-1].to = -1;
}

///////////////////////////////////////////////////////////////////////////////////////////
list::~list()
{
	delete[] itmlist;
}

///////////////////////////////////////////////////////////////////////////////////////////
int list::initialize(double x)
{
	ptr = 0;
	start = 0;
	stop = len-1;
	int i;
	for(i=0;i<len;i++)
	{
		itmlist[i].value = x;
		itmlist[i].from = i-1;
		itmlist[i].to = i+1;
	}
	itmlist[len-1].to = -1;
	
	/*	cout<<'\n'<<ptr<<'\n'<<"index"<<'\t'<<"value"<<'\t'<<"from"<<'\t'<<"to"<<'\t'<<"start"<<'\t'<<"stop"<<endl;
	for(i=0;i<len;i++)
	{
	cout<<i<<'\t'<<itmlist[i].value<<'\t'<<itmlist[i].from<<'\t'<<itmlist[i].to<<'\t'<<start<<'\t'<<stop<<endl;
	}
	*/	
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
int list::insert(double data)
{
	itmlist[ptr].value = data;
	this->resort();
	ptr = (ptr+1)%len;
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
int list::resort(void)
{
	int i,temp;
	
	// Moving forward in the list
	while(itmlist[ptr].to> -1 && itmlist[ptr].value>itmlist[itmlist[ptr].to].value)
	{
		temp = itmlist[ptr].to;
		if(itmlist[temp].to > -1)
		{
			itmlist[itmlist[temp].to].from = ptr;
		}
		itmlist[ptr].to = itmlist[temp].to;
		itmlist[temp].from = itmlist[ptr].from;
		if(itmlist[ptr].from > -1)
		{
			itmlist[itmlist[ptr].from].to = temp;
		}
		itmlist[ptr].from = temp;
		itmlist[temp].to = ptr;
	}
	// Moving forward in the list
	while(itmlist[ptr].from> -1 && itmlist[ptr].value<itmlist[itmlist[ptr].from].value)
	{
		temp = itmlist[ptr].from;
		if(itmlist[temp].from > -1)
		{
			itmlist[itmlist[temp].from].to = ptr;
		}
		itmlist[ptr].from = itmlist[temp].from;
		itmlist[temp].to = itmlist[ptr].to;
		if(itmlist[ptr].to > -1)
		{
			itmlist[itmlist[ptr].to].from = temp;
		}
		itmlist[ptr].to = temp;
		itmlist[temp].from = ptr;
	}
	
	for(i=0;i<len;i++)
	{
		if(itmlist[i].from==-1)
		{
			start = i;
			break;
		}
	}
	
	for(i=0;i<len;i++)
	{
		if(itmlist[i].to==-1)
		{
			stop = i;
			break;
		}
	}
	
	
	/*	cout<<'\n'<<ptr<<'\n'<<"index"<<'\t'<<"value"<<'\t'<<"from"<<'\t'<<"to"<<'\t'<<"start"<<'\t'<<"stop"<<endl;
	for(i=0;i<len;i++)
	{
	cout<<i<<'\t'<<itmlist[i].value<<'\t'<<itmlist[i].from<<'\t'<<itmlist[i].to<<'\t'<<start<<'\t'<<stop<<endl;
	}
	*/	
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////
double list::mean(int a, int b) // neglects the first 'a' and last 'b' samples 
{
	int pntr = start;
	int cntr = 0;
	double sum = 0;
	while(cntr<len-b){
		if(cntr>=a){
			sum += itmlist[pntr].value;
		}
		cntr++;
		pntr = itmlist[pntr].to;
	}

	return sum/(cntr-a);
}

///////////////////////////////////////////////////////////////////////////////////////////
double list::wmedian(const double *weight, int m, int a, int b) // neglects the first 'a' and last 'b' samples 
{
	int pntr = start;
	int cntr = 0;
	double sum = 0;
	while(cntr<len-b){
		if(cntr>=a){
			sum += weight[cntr]*itmlist[pntr].value;
		}
		cntr++;
		pntr = itmlist[pntr].to;
	}

	return sum;
}

///////////////////////////////////////////////////////////////////////////////////////////
double list::medianOdd(void)
{
	int pntr = start;
	int cntr = 0;
	while(cntr<len/2)
	{
		pntr = itmlist[pntr].to;
		cntr++;
	}
	return itmlist[pntr].value;
}

///////////////////////////////////////////////////////////////////////////////////////////
double list::medianEven(void)
{
	int pntr = start;
	int cntr = 0;
	while(cntr<len/2)
	{
		pntr = itmlist[pntr].to;
		cntr++;
	}
	return (itmlist[pntr].value + itmlist[itmlist[pntr].from].value)/2.;
}

///////////////////////////////////////////////////////////////////////////////////////////
