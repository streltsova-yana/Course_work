#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <vector>

#include "global_search.h"

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

void segment_split(double (*fcnPtr)(double), double a, double b, int kmax, double precision, double& xmin, double& ymin)
{
    if (a > b)
        throw "Incorrect bounds: a must be less than b";

    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double segment = static_cast<double>(b - a) / size;
    double local_a = a + segment * rank;
    double local_b = a + segment * (rank + 1);

    struct {
        double ymin = std::numeric_limits<double>::max();
        double xmin = std::numeric_limits<double>::max();
    }local_min, global_min;

    sequential_global_search(fcnPtr, local_a, local_b, kmax, precision, local_min.xmin, local_min.ymin);

    MPI_Reduce(&local_min, &global_min, 1, MPI_2DOUBLE_PRECISION, MPI_MINLOC, 0, MPI_COMM_WORLD);

    xmin = global_min.xmin;
    ymin = global_min.ymin;
}

void parallel_operations(double (*fcnPtr)(double), double a, double b, int kmax, double precision, double& xmin, double& ymin)
{
    if (a > b)
        throw "Incorrect bounds: a must be less than b";

    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<std::pair<double, double>> x = { {a, fcnPtr(a)}, {b, fcnPtr(b)} };
    int k = 2;
    bool prec = 1;
    while ((k < kmax) && prec) {
        double local_M = 0, global_M;
        for (int i1 = rank, i2 = rank + 1; i2 < x.size(); i1 += size, i2 += size) {
            local_M = fmax(local_M, abs((x[i2].second - x[i1].second) / (x[i2].first - x[i1].first)));
        }
        MPI_Allreduce(&local_M, &global_M, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double r = 2;
        double m = r * global_M;
        if (global_M == 0)
            m = 1;

        struct {
            double val = 0;
            int pos = 0;
        }local_R, global_R;
        for (int i1 = rank, i2 = rank + 1; i2 < x.size(); i1 += size, i2 += size) {
            double Ri = m * (x[i2].first - x[i1].first) + (x[i2].second - x[i1].second) * (x[i2].second - x[i1].second) /
                (m * (x[i2].first - x[i1].first)) - 2 * (x[i2].second - x[i1].second);
            if (Ri > local_R.val) {
                local_R.val = Ri;
                local_R.pos = i1;
            }
        }
        MPI_Allreduce(&local_R, &global_R, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

        double new_x = 0.5 * (x[global_R.pos + 1].first + x[global_R.pos].first) - (x[global_R.pos + 1].second - x[global_R.pos].second) / (2 * m);
        x.emplace(x.begin() + global_R.pos + 1, std::pair<double, double>(new_x, fcnPtr(new_x)));
        k++;
        prec = x[global_R.pos + 1].first - x[global_R.pos].first <= precision ? 0 : 1;
    }
    xmin = x[0].first;
    ymin = x[0].second;
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

    int size, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<std::pair<double, double>> x = { {a, fcnPtr(a)}, {b, fcnPtr(b)} };
    int k = 2;
    int prec = 1;
    while ((k < kmax) && prec) {

        double local_M = 0, global_M = 0;
        for (int i1 = rank, i2 = rank + 1; i2 < x.size(); i1 += size, i2 += size) {
            local_M = fmax(local_M, abs((x[i2].second - x[i1].second) / (x[i2].first - x[i1].first)));
        }
        MPI_Allreduce(&local_M, &global_M, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        double r = 2;
        double m = r * global_M;
        if (global_M == 0)
            m = 1;

        struct R {
            double val = 0;
            int pos = 0;
        };
        std::vector<R> local_Ri;
        for (int i1 = rank, i2 = rank + 1; i2 < x.size(); i1 += size, i2 += size) {
            local_Ri.push_back(R{ m * (x[i2].first - x[i1].first) + (x[i2].second - x[i1].second) * (x[i2].second - x[i1].second) /
                (m * (x[i2].first - x[i1].first)) - 2 * (x[i2].second - x[i1].second), i1});
        }
        
        std::vector<R> global_Ri(x.size() - 1);
        int part = (x.size() - 1) / size;
        std::vector<int> scounts(size, part);
        for (int i = 0; i < (x.size() - 1) % size; i ++)
            scounts[i]++;
        std::vector<int> displs(size, 0);
        for (int i = 1; i < size; i++)
            displs[i] = displs[i - 1] + scounts[i - 1];
        MPI_Allgatherv(local_Ri.data(), local_Ri.size(), MPI_DOUBLE_INT, global_Ri.data(), scounts.data(), displs.data(), MPI_DOUBLE_INT, MPI_COMM_WORLD);
        int n = global_Ri.size() < size ? global_Ri.size() : size;
        std::sort(global_Ri.begin(), global_Ri.end(), [](R a, R b)->bool {return a.val > b.val; });
        
        struct Point {
            double x;
            double y;
        } local_point;
        int local_prec = prec;
        for (int i = rank; i < n; i += size) {
            local_point.x = 0.5 * (x[global_Ri[i].pos + 1].first + x[global_Ri[i].pos].first) - (x[global_Ri[i].pos + 1].second - x[global_Ri[i].pos].second) / (2 * m);
            local_point.y = fcnPtr(local_point.x);
            if (local_prec)
                local_prec = x[global_Ri[i].pos + 1].first - x[global_Ri[i].pos].first <= precision ? 0 : 1;
        }
        MPI_Allreduce(&local_prec, &prec, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);
        
        std::fill(scounts.begin(), scounts.end(), n / size);
        for (int i = 0; i < n % size; i++)
            scounts[i]++;
        for (int i = 1; i < size; i++)
            displs[i] = displs[i - 1] + scounts[i - 1];
        std::vector<Point> global_point(n);
        MPI_Allgatherv(&local_point, scounts[rank], MPI_2DOUBLE_PRECISION, global_point.data(), scounts.data(), displs.data(), MPI_2DOUBLE_PRECISION, MPI_COMM_WORLD);

        for (int i = 0; i < n; i++)
            x.push_back(std::pair<double, double>(global_point[i].x, global_point[i].y));
        std::sort(x.begin(), x.end(), [](std::pair<double, double> a, std::pair<double, double> b)->bool {return a.first < b.first; });
        k++;
    }
    xmin = x[0].first;
    ymin = x[0].second;
    for (const auto& i : x) {
        if (ymin > i.second) {
            ymin = i.second;
            xmin = i.first;
        }
    }
}