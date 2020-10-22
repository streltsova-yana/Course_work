#ifndef GLOBAL_SEARCH_H
#define GLOBAL_SEARCH_H

double func1(double x);

double func2(double x);

void global_search(double (*fcnPtr)(double), double a, double b, int kmax, double precision, double &xmin, double &ymin);

#endif