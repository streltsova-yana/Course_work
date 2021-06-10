#include "tbb/tbb.h"
#define TBB_PREVIEW_CONCURRENT_ORDERED_CONTAINERS 1
#include "tbb/concurrent_map.h"

#include <map>
#include <cmath>
#include <iterator>

#include "global_search.h"

using namespace tbb;

void sequential_global_search(double (*fcnPtr)(double), double a, double b, int kmax, double precision, double& xmin, double& ymin)
{
    if (a > b)
        throw "Incorrect bounds: a must be less than b";
    std::map<double, double> x = { {a, fcnPtr(a)}, {b, fcnPtr(b)} };
    int k = 2;
    bool prec = 1;
    while ((k < kmax) && prec) {
        double M = 0;
        for (auto it1 = x.begin(), it2 = ++x.begin(); it2 != x.end(); it1++, it2++) {
            M = fmax(M, abs((it2->second - it1->second) / (it2->first - it1->first)));
        }
        double r = 2;
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
        prec = i2 - i1 <= precision ? 0 : 1;
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

void parallel_global_search(double (*fcnPtr)(double), double a, double b, int kmax, double precision, double& xmin, double& ymin)
{
    if (a > b)
        throw "Incorrect bounds: a must be less than b";
    concurrent_map<double, double> x = { {a, fcnPtr(a)}, {b, fcnPtr(b)} };
    int k = 2;
    bool prec = 1;
    while ((k < kmax) && prec) {
        double M = 0;
        for (auto it1 = x.begin(), it2 = ++x.begin(); it2 != x.end(); it1++, it2++) {
            M = fmax(M, abs((it2->second - it1->second) / (it2->first - it1->first)));
        }
        double r = 2;
        double m = r * M;
        if (M == 0)
            m = 1;
        concurrent_map<double, std::pair<double, double>, std::greater<double>> R;
        int grainSize = x.size() < 4 ? 1 : x.size() / 4;
        parallel_for(blocked_range<size_t>(0, x.size(), grainSize),
            [&](const blocked_range<size_t>& r) {
                for (auto it1 = x.begin(), it2 = ++x.begin(); it2 != x.end(); it1++, it2++) {
                    double Ri = m * (it2->first - it1->first) + (it2->second - it1->second) * (it2->second - it1->second) /
                        (m * (it2->first - it1->first)) - 2 * (it2->second - it1->second);
                    R[Ri] = std::pair<double, double>(it1->first, it2->first);
                }
            });

        int size = R.size() < 4 ? R.size() : 4;
        parallel_for(blocked_range<size_t>(0, size, 1),
            [&](const blocked_range<size_t>& r) {
                auto begin = R.begin();
                auto end = R.begin(); 
                std::advance(end, size);
                for (auto it = begin; it != end; it++) {
                    prec = it->second.second - it->second.first <= precision ? 0 : 1;
                    double new_x = 0.5 * (it->second.second + it->second.first) - (x[it->second.second] - x[it->second.first]) / (2 * m);
                    x[new_x] =  fcnPtr(new_x);
                }
            });

        k++;
    }
    xmin = a;
    ymin = x[a];
    for (const auto& i : x) {
        if (ymin > i.second) {
            ymin = i.second;
            xmin = i.first;
        }
    }
};
