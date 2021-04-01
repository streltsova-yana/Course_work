#include <iostream>

#include "Segment_split.h"
#include "hansen_functions.h"

using namespace std;
using namespace tbb;

int main()
{
    int kmax = 1000;
    double precision = 0.001;
    double xmin, ymin;
    double a, b;
    int grainSize;

    for (size_t i = 0; i < 20; i++) {
        xmin = numeric_limits<double>::max();
        ymin = numeric_limits<double>::max();
        a = intervals[i][0];
        b = intervals[i][1];
        grainSize = (b - a) / 4;
        if (grainSize <= 0)
            grainSize = 1;
        task_scheduler_init(4);
        try {
            Segment_split s(pfn[i], kmax, precision);
            parallel_reduce(blocked_range<double>(a, b, grainSize), s);
            s.result(xmin, ymin);
            cout << "Global_min_hfunc " << i + 1 << " - " << res[i][0]
                 << ";     Result: (" << xmin << ", " << ymin << ")" <<endl;
        }
        catch (const char* message) {
            cout << message;
        }

    }
}