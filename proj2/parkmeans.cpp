#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

const int MAX_ITER = 100;
const int n_means = 4;

int file_size(ifstream& file){
    file.seekg(0, ios::end);
    int n = file.tellg();
    file.seekg(0, ios::beg);
    return n;
}

double distance(double *v1, double *v2, int n)
{
    double d = 0;
    for (int i = 0; i < n; i++) {
        d += (v1[i] - v2[i]) * (v1[i] - v2[i]);
    }
    return sqrt(d);
}

int find_nearest_cluster(double *value, double **clusters, int n, int m)
{
    int cluster_index = -1;
    double min_distance = INFINITY;
    for (int i = 0; i < n; i++) {
        double d = distance(value, clusters[i], m);
        if (d < min_distance) {
            min_distance = d;
            cluster_index = i;
        }
    }
    return cluster_index;
}

void compute_clusters(double **values, int n, int m, double **clusters, int *assignments)
{
    for (int i = 0; i < n_means; i++) {
        int idx = rand() % n;
        for (int j = 0; j < m; j++) {
            clusters[i][j] = values[idx][j];
        }
    }

    bool changed = true;
    int iter = 0;
    while (changed && iter < MAX_ITER) {
        changed = false;
        for (int i = 0; i < n; i++) {
            int nearest_cluster = find_nearest_cluster(values[i], clusters, n_means, m);
            if (nearest_cluster != assignments[i]) {
                changed = true;
                assignments[i] = nearest_cluster;
            }
        }

        for (int i = 0; i < n_means; i++) {
            double *sum = new double[m];
            int count = 0;
            for (int j = 0; j < n; j++) {
                if (assignments[j] == i) {
                    count++;
                    for (int k = 0; k < m; k++) {
                        sum[k] += values[j][k];
                    }
                }
            }
            if (count > 0) {
                for (int k = 0; k < m; k++) {
                    clusters[i][k] = sum[k] / count;
                }
            }
            delete[] sum;
        }

        iter++;
    }
}

int main(int argc, char *argv[])
{

    int rank;
    int nprocs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    //number of processes must be 4 - 32
    if (nprocs < 4 || nprocs > 32){
        std::cerr << "Wrong ammount of processes" << std::endl;
        MPI_Finalize();
        return 1;
    }

    //opening binary file
    std::ifstream input("numbers", std::ios::binary);
    
    if (!input){
        std::cerr << "Failed to open file" << std::endl;
        MPI_Finalize();
        return 1;
    }
    int n = file_size(input); 
    printf("%d", n);

    int m = 1;   // number of dimensions

    // Allocate memory for the values and assignments
    double **values = new double *[n];
    int *assignments = new int[n];

    for (int i = 0; i < n; i++) {
        values[i] = new double[m];
        assignments[i] = -1;
    }

    // Generate random values
    srand(12345);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            values[i][j] = (double) rand() / RAND_MAX;
        }
    }

    // Allocate memory for the clusters
    double **clusters = new double *[n_means];
    for (int i = 0; i < n_means; i++) {
        clusters[i] = new double[m];
    }

    // Compute the clusters using K-Means
    compute_clusters(values, n, m, clusters, assignments);

    // Print the results
    if (rank == 0) {
        printf("Clusters:\n");
        for (int i = 0; i < n_means; i++) {
            for (int j = 0; j < m; j++) {
                printf("%f ", clusters[i][j]);
            }
            printf("\n");
        }
    }

    // Deallocate memory
    for (int i = 0; i < n; i++) {
        delete[] values[i];
    }
    delete[] values;
    delete[] assignments;
    for (int i = 0; i < n_means; i++) {
        delete[] clusters[i];
    }
    delete[] clusters;

    // Finalize MPI
    MPI_Finalize();

    return 0;
}



