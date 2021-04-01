#ifndef __SEGMENT_SPLIT_H__
#define __SEGMENT_SPLIT_H__
#define NOMINMAX

#include <tbb/tbb.h>

#include <limits>

#include "global_search.h"

using namespace std;
using namespace tbb;

class Segment_split
{
    double (*fcnPtr)(double);
    int kmax;
    double precision;
    double xmin, ymin;
public:
    explicit Segment_split(double (*_fcnPtr)(double), int _kmax, double _precision) :
        fcnPtr(_fcnPtr), kmax(_kmax), precision(_precision), xmin(numeric_limits<double>::max()), ymin(numeric_limits<double>::max()) {};

    Segment_split(const Segment_split& tmp, split) :
        fcnPtr(tmp.fcnPtr), kmax(tmp.kmax), precision(tmp.precision), xmin(numeric_limits<double>::max()), ymin(numeric_limits<double>::max()) {};

    void operator()(const blocked_range<double>& r) {
        double begin = r.begin(), end = r.end();
        global_search(fcnPtr, begin, end, kmax, precision, xmin, ymin);
    };

    void join(const Segment_split& tmp) {
        if (ymin > tmp.ymin) {
            ymin = tmp.ymin;
            xmin = tmp.xmin;
        }
    };

    void result(double& _xmin, double& _ymin) {
        _xmin = xmin;
        _ymin = ymin;
    };
};

#endif
