#include <iostream>
#include <limits>
#include "global_search.h"

using namespace std;

int main()
{
	int kmax = 1000;
	double precision = 0.1;
	double xmin = numeric_limits<double>::max(), ymin = numeric_limits<double>::max();
	double a = -73, b = 45.5;
	try {
		global_search(func1, a, b, kmax, precision, xmin, ymin);
		cout << "Global_min_func1: (" << xmin << ", " << ymin << ")" << endl;
	}
	catch (const char* message) {
		cout << message;
	}
	a = 1.2;
	b = 3.14;
	try {
		global_search(func2, a, b, kmax, precision, xmin, ymin);
		cout << "Global_min_func2: (" << xmin << ", " << ymin << ")" << endl;
	}
	catch (const char* message) {
		cout << message;
	}
	a = 1.2;
	b = 3.14;
	try {
		global_search(func2, b, a, kmax, precision, xmin, ymin);
		cout << "Global_min_func2: (" << xmin << ", " << ymin << ")" << endl;
	}
	catch (const char* message) {
		cout << message;
	}
}