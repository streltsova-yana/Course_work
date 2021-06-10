#ifndef __GLOBAL_SEARCH_H__
#define __GLOBAL_SEARCH_H__

void sequential_global_search(double (*fcnPtr)(double), double a, double b, int kmax, double precision, double& xmin, double& ymin);

void segment_split(double (*fcnPtr)(double), double a, double b, int kmax, double precision, double& xmin, double& ymin);

void parallel_operations(double (*fcnPtr)(double), double a, double b, int kmax, double precision, double& xmin, double& ymin);

void parallel_global_search(double (*fcnPtr)(double), double a, double b, int kmax, double precision, double& xmin, double& ymin);

#endif
