#ifndef __GLOBAL_SEARCH_H__
#define __GLOBAL_SEARCH_H__

void global_search(double (*fcnPtr)(double), double a, double b, int kmax, double precision, double& xmin, double& ymin);

#endif