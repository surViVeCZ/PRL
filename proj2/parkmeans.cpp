#include <mpi.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>

const int MAX_ITER = 100;
const int N_CLUSTERS = 4;
const int MAX_VALUES = 32;

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
    for (int i = 0; i < N_CLUSTERS; i++) {
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
            int nearest_cluster = find_nearest_cluster(values[i], clusters, N_CLUSTERS, m);
            if (nearest_cluster != assignments[i]) {
                changed = true;
                assignments[i] = nearest_cluster;
            }
        }

        for (int i = 0; i < N_CLUSTERS; i++) {
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
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get the rank of the current process
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the number of processes
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Set the number of values and dimensions
    int n = 16;  // number of values
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
    double **clusters = new double *[N_CLUSTERS];
    for (int i = 0; i < N_CLUSTERS; i++) {
        clusters[i] = new double[m];
    }

    // Compute the clusters using K-Means
    compute_clusters(values, n, m, clusters, assignments);

    // Print the results
    if (rank == 0) {
        printf("Clusters:\n");
        for (int i = 0; i < N_CLUSTERS; i++) {
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
    for (int i = 0; i < N_CLUSTERS; i++) {
        delete[] clusters[i];
    }
    delete[] clusters;

    // Finalize MPI
    MPI_Finalize();

    return 0;
}



