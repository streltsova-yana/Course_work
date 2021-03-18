#include <map>
#include <cmath>

#include "global_search.h"

void global_search(double (*fcnPtr)(double), double a, double b, int kmax, double precision, double& xmin, double& ymin)
{
	if (a > b)
		throw "Incorrect bounds: a must be less than b";
	std::map<double, double> x = { {a, fcnPtr(a)}, {b, fcnPtr(b)} };
	int k = 2;
	bool prec = 1;
	while ((k < kmax) && prec) {
		double M = 0;
		for (auto it1 = x.begin(), it2 = ++x.begin(); it2 != x.end(); it1++, it2++)
			M = fmax(M, abs((it2->second - it1->second) / (it2->first - it1->first)));
		double r = 100;
		double m = r * M;
		if (M == 0)
			m = 1;
		double R = 0, i1 = 0, i2 = 0;
		for (auto it1 = x.begin(), it2 = ++x.begin(); it2 != x.end(); it1++, it2++) {
			double Ri = m * (it2->first - it1->first) + (it2->second - it1->second) * (it2->second - it1->second) /
				(m * (it2->first - it1->first)) - 2 * (it2->second - it1->second);
			if (Ri > R) {
				R = Ri;
				i1 = it1->first;
				i2 = it2->first;
			}
		}
		double new_x = 0.5 * (i2 + i1) - (x[i2] - x[i1]) / (2 * m);
		x[new_x] = fcnPtr(new_x);
		k++;
		if (i2 - i1 <= precision)
			prec = 0;
	}
	xmin = a;
	ymin = x[a];
	for (const auto& i : x) {
		if (ymin > i.second) {
			ymin = i.second;
			xmin = i.first;
		}
	}
}