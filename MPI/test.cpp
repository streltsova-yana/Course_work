#include <mpi.h>

#include <iostream>
#include <ctime>

#include "hansen_functions.h"
#include "global_search.h"

using namespace std;

int main(int argc, char* argv[])
{
    int kmax = 1000;
    double precision = 0.01;
    double xmin_par_split, ymin_par_split;
    double xmin_par_oper, ymin_par_oper;
    double xmin_par, ymin_par;
    double xmin_seq, ymin_seq;
    double a, b;
    int grainSize;

    MPI_Init(&argc, &argv);
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (size_t i = 0; i < 20; i++) {
        xmin_par_split = numeric_limits<double>::max();
        ymin_par_split = numeric_limits<double>::max();
        xmin_par_oper = numeric_limits<double>::max();
        ymin_par_oper = numeric_limits<double>::max();
        xmin_par = numeric_limits<double>::max();
        ymin_par = numeric_limits<double>::max();
        xmin_seq = numeric_limits<double>::max();
        ymin_seq = numeric_limits<double>::max();
        a = intervals[i][0];
        b = intervals[i][1];

        try {
            double t1_par_split = MPI_Wtime();
            segment_split(pfn[i], a, b, kmax, precision, xmin_par_split, ymin_par_split);
            double t2_par_split = MPI_Wtime();

            double t1_par_oper = MPI_Wtime();
            parallel_operations(pfn[i], a, b, kmax, precision, xmin_par_oper, ymin_par_oper);
            double t2_par_oper = MPI_Wtime();

            double t1_par = MPI_Wtime();
            parallel_global_search(pfn[i], a, b, kmax, precision, xmin_par, ymin_par);
            double t2_par = MPI_Wtime();


            if (rank == 0) {
                double t1_seq = MPI_Wtime();
                sequential_global_search(pfn[i], a, b, kmax, precision, xmin_seq, ymin_seq);
                double t2_seq = MPI_Wtime();

                cout << "Global_min_hfunc " << i + 1 << " - " << res[i][0] << endl;
                cout << "Result_parallel_split: (" << xmin_par_split << ", " << ymin_par_split << ")     Time_parallel: " << t2_par_split - t1_par_split << endl;
                cout << "Result_parallel_operations: (" << xmin_par_oper << ", " << ymin_par_oper << ")     Time_parallel: " << t2_par_oper - t1_par_oper << endl;
                cout << "Result_parallel: (" << xmin_par << ", " << ymin_par << ")     Time_parallel: " << t2_par - t1_par << endl;
                cout << "Result_sequential: (" << xmin_seq << ", " << ymin_seq << ")   Time_sequential: " << t2_seq - t1_seq << endl << endl;
            }
        }
        catch (const char* message) {
            cout << message;
        }

    }
    MPI_Finalize();
}
