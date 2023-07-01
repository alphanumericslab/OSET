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

#define SWAP(a,b) temp = (a);(a) = (b);(b) = temp; 

class item
{
public:
	double value;
	int from;
	int to;
};

class list
{
private:
	item* itmlist;
	int ptr;
	int len;
	int start,stop;

	int resort(void);

public:
	list(int N,double init=0);
	~list(void);
	
	int initialize(double x);
	int insert(double data);
	
	double mean(int alpha);
	double wmedian(const double *weight, int m, int alpha);
	double medianOdd(void);
	double medianEven(void);
};