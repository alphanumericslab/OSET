#ifndef _LIST_H_
#define _LIST_H_
// A sorted list class used for trimmed mean and median filtering

/*
% Revision History:
%   2008: First release
%   2023: Renamed from deprecated version List.h
%
% Reza Sameni, 2022-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

*/

#define SWAP(a,b) temp = (a);(a) = (b);(b) = temp; 

class item
{
public:
	double value;
	int value_index;
	int from;
	int to;
};

class list
{
private:
	item* itmlist;
	int write_pointer;
	int list_len;
	int head_pointer;
	int tail_pointer;

	int resort(void);

public:
	list(int N, double init = 0, int init_index = 0);
	~list(void);
	
	int initialize(double x, int init_index = 0);
	int insert(double data);
	int insert(double data, int data_index);
	
	double mean(int alpha);
	double wmedian(const double *weight, int m, int alpha);
	double median_odd(void);
	double median_even(void);
	double max(void);
	int max_index(void);
	double min(void);
	int min_index(void);
};
#endif /*_LIST_H_*/