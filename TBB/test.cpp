#include <iostream>
#include <ctime>

#include "Segment_split.h"
#include "hansen_functions.h"

using namespace std;
using namespace tbb;

int main()
{
    int kmax = 1000;
    double precision = 0.01;
    double xmin_par_split, ymin_par_split;
    double xmin_par, ymin_par;
    double xmin_seq, ymin_seq;
    double a, b;
    int grainSize;
    task_scheduler_init init(4);

    for (size_t i = 0; i < 20; i++) {
        xmin_par_split = numeric_limits<double>::max();
        ymin_par_split = numeric_limits<double>::max();
        xmin_par = numeric_limits<double>::max();
        ymin_par = numeric_limits<double>::max();
        xmin_seq = numeric_limits<double>::max();
        ymin_seq = numeric_limits<double>::max();
        a = intervals[i][0];
        b = intervals[i][1];
        grainSize = (b - a) / 4;
        if (grainSize <= 0)
            grainSize = 1;

        try {
            tick_count t1_par_split = tick_count::now();
            Segment_split s(pfn[i], kmax, precision);
            parallel_reduce(blocked_range<double>(a, b, grainSize), s);
            tick_count t2_par_split = tick_count::now();
            s.result(xmin_par_split, ymin_par_split);

            tick_count t1_par = tick_count::now();
            parallel_global_search(pfn[i], a, b, kmax, precision, xmin_par, ymin_par);
            tick_count t2_par = tick_count::now();

            tick_count t1_seq = tick_count::now();
            sequential_global_search(pfn[i], a, b, kmax, precision, xmin_seq, ymin_seq);
            tick_count t2_seq = tick_count::now();

            cout << "Global_min_hfunc " << i + 1 << " - " << res[i][0] << endl;
            cout << "Result_parallel_split: (" << xmin_par_split << ", " << ymin_par_split << ")     Time_parallel: " << (t2_par_split - t1_par_split).seconds() << endl;
            cout << "Result_parallel: (" << xmin_par << ", " << ymin_par << ")     Time_parallel: " << (t2_par - t1_par).seconds() << endl;
            cout << "Result_sequential: (" << xmin_seq << ", " << ymin_seq << ")   Time_sequential: " << (t2_seq - t1_seq).seconds() << endl << endl;
        }
        catch (const char* message) {
            cout << message;
        }

    }
}
