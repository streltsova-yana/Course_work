#include <iostream>
#include <limits>
#include "global_search.h"

using namespace std;

int main()
{
	int kmax = 100;
	double precision = 1e-5;
	double xmin = numeric_limits<double>::max(), ymin = numeric_limits<double>::max();
	double a = -1, b = 2;
	cout << endl;
	try {
		global_search(squaring, a, b, kmax, precision, xmin, ymin);
		cout << "Global_min_squaring: (" << xmin << ", " << ymin
			 << ") on the segment [" << a << ", " << b << "]" << endl;
	}
	catch (const char* message) {
		cout << message;
	}
	a = 1.2;
	b = 3.14;
	try {
		global_search(sinus, a, b, kmax, precision, xmin, ymin);
		cout << "Global_min_sinus: (" << xmin << ", " << ymin
			 << ")" << ") on the segment [" << a << ", " << b << "]" << endl;
	}
	catch (const char* message) {
		cout << message;
	}
	a = 1.2;
	b = 3.14;
	try {
		global_search(sinus, b, a, kmax, precision, xmin, ymin);
		cout << "Global_min_sinus: (" << xmin << ", " << ymin
			 << ")" << ") on the segment [" << a << ", " << b << "]" << endl;
	}
	catch (const char* message) {
		cout << message << endl;
	}
}