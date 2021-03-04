#ifndef GLOBAL_SEARCH_H
#define GLOBAL_SEARCH_H

double squaring(double x);

double sinus(double x);

void global_search(double (*fcnPtr)(double), double a, double b, int kmax, double precision, double &xmin, double &ymin);

#endif