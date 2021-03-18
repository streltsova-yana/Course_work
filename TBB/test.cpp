#include <iostream>

#include "Segment_split.h"
#include "hansen_functions.h"

using namespace std;
using namespace tbb;

int main()
{
	int kmax = 1000;
	double precision = 0.1;
	double xmin = numeric_limits<double>::max(), ymin = numeric_limits<double>::max();
	double a = -1.5, b = 11;//intervals[0][1], b = intervals[0][2];
	int grainSize = (b - a) / 4;
	task_scheduler_init(4);
	try {
		Segment_split s(hfunc1, kmax, precision);
		parallel_reduce(blocked_range<double>(a, b, grainSize), s);
		s.result(xmin, ymin);
		cout << "Global_min_hfunc1: (" << xmin << ", " << ymin << ")" << endl;
	}
	catch (const char* message) {
		cout << message;
	}
}