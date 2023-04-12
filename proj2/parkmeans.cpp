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
#include <iomanip>

using namespace std;
const int n_means = 4;

double distance(unsigned char value, double mean) {
    return abs(value - mean);
}

int find_nearest_cluster(unsigned char value, const vector<double>& means) {
    int cluster_index = -1;
    double min_distance = INFINITY;
    for (int i = 0; i < n_means; i++) {
        double d = distance(value, means[i]);
        if (d < min_distance) {
            min_distance = d;
            cluster_index = i;
        }
    }
    return cluster_index;
}

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


void print_cluster(int i, const vector<double>& means, const vector<vector<unsigned char> >& clusters)
{
    std::cout << "Cluster " << i+1 << ": ";
    std::cout << std::setprecision(1) << std::fixed << "[" << means[i] << "]" << " ";
    for (int j = 0; j < clusters[i].size(); j++) {
        //for last dont print comma
        if( j == clusters[i].size()-1)
            std::cout << static_cast<int>(clusters[i][j]);
        else
            std::cout << static_cast<int>(clusters[i][j]) << ", ";
    }
    std::cout << std::endl;
}

//calculating new clusters based on distances, returning assignments
std::vector<int> new_clusters(int process_no, int n, const vector<unsigned char>& numbers_file, const vector<int>& assignments)
{
 //create new clusters from numbers_file based on assignments
    std::vector<std::vector<unsigned char>> clusters(n_means);
    for (int i = 0; i < n; i++) {
        clusters[assignments[i]].push_back(numbers_file[i]);
    }
    //calculate new means
    std::vector<double> means(n_means);
    for (int i = 0; i < n_means; i++) {
        double sum = 0;
        for (int j = 0; j < clusters[i].size(); j++) {
            sum += clusters[i][j];
        }
        if (clusters[i].size() > 0) {
            means[i] = sum / clusters[i].size();
        } else {
            means[i] = 0;
        }
    }
    // print clusters
    if (process_no == 0) {
        for (int i = 0; i < n_means; i++) {
            print_cluster(i, means, clusters);
        }
    }
    //calculate new assignments
    std::vector<int> new_assignments(n);
    for (int i = 0; i < n; i++) {
        new_assignments[i] = find_nearest_cluster(numbers_file[i], means);
    }
    if (process_no == 0) {
        // std::cout << "Assignments: ";
        // for (int i = 0; i < n; i++) {
        //     std::cout << assignments[i] << " ";
        // }
        std::cout << std::endl;
    }
    return new_assignments;
}


#include <mpi.h>
#include <vector>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[]) {
    int process_no;
    int nprocs;
    std::vector<double> means(n_means);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_no);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    std::ifstream input("numbers", std::ios::binary);
    int n = 0;
    if (process_no == 0) {
        n = file_size(input);

        if (nprocs > n) {
            std::cerr << "Too many processes" << std::endl;
            MPI_Finalize();
            return 1;
        }
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (n > nprocs) {
        input.seekg(n-nprocs, std::ios::beg);
        n = nprocs;
    }

    std::vector<unsigned char> numbers_file(n);
    std::vector<double> all_means(n_means * nprocs);

    if (process_no == 0) {
        input.read(reinterpret_cast<char*>(&numbers_file[0]), n);
        input.close();
        printf("Initialized centers: ");
        for (int i = 0; i < n_means; i++) {
            means[i] = numbers_file[i];
            printf("%d ", int(means[i]));
        }
        std::cout << std::endl;
    }

    MPI_Scatter(&numbers_file[0], 1, MPI_UNSIGNED_CHAR, &numbers_file[0], 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(&means[0], n_means, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //declare assignments
    std::vector<int> assignments(n);
    for (int i = 0; i < n; i++) {
        assignments[i] = find_nearest_cluster(numbers_file[i], means);
    }


    //calculate new clusters based on distances, returning assignments, do it parallel using MPI gatherall and reduce
    while(true) {
        std::vector<int> new_assignments = new_clusters(process_no, n, numbers_file, assignments);
        if (new_assignments == assignments) {
            break;
        }
        assignments = new_assignments;
        //compute new means and send it to all processes
        for (int i = 0; i < n_means; i++) {
            double sum = 0;
            int count = 0;
            for (int j = 0; j < n; j++) {
                if (assignments[j] == i) {
                    sum += numbers_file[j];
                    count++;
                }
            }
            if (count > 0) {
                means[i] = sum / count;
            } else {
                means[i] = 0;
            }
        }
        MPI_Allgather(&means[0], n_means, MPI_DOUBLE, &all_means[0], n_means, MPI_DOUBLE, MPI_COMM_WORLD);
    }
    //print new means of all processes
    if (process_no == 0) {
        printf("Final centers: ");
        for (int i = 0; i < n_means; i++) {
            printf("%d ", int(means[i]));
        }
        std::cout << std::endl;
    }



    MPI_Finalize();
    return 0;
}




