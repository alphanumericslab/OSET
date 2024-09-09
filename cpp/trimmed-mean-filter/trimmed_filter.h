#ifndef _TRIMMED_FILTER_H_
#define _TRIMMED_FILTER_H_

int trimmed_filter(const double *arr, double *output, int n, const char* type, int w, int alpha = 0, const double* h = 0);
int trimmed_filter_max_min_indexes(const double *arr, int *output, int n, const char* type, int w);

#endif /*_TRIMMED_FILTER_H_*/