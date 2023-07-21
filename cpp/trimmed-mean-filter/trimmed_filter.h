#ifndef _TRIMMED_FILTER_H_
#define _TRIMMED_FILTER_H_

int trimmed_filter(const double *arr, double *output, int n, const char* type, int w, int alpha = 0, const double* h = 0);

#endif /*_TRIMMED_FILTER_H_*/