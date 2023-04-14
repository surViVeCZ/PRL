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

int find_nearest_cluster(int value, const vector<double>& means) {
    int cluster_index = -1;
    double min_distance = INFINITY;
    for (int i = 0; i < n_means; i++) {
        double d = abs(value - means[i]);
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


void print_cluster(double mean, const vector<int >& cluster)
{
    std::cout << std::setprecision(1) << std::fixed << "[" << mean << "]" << " ";
    for (int j = 0; j < cluster.size(); j++) {
        //for last dont print comma
        if( j == cluster.size()-1)
            std::cout << cluster[j];
        else
            std::cout << cluster[j] << ", ";
    }
    std::cout << std::endl;
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

    std::vector<int> numbers(nprocs);
    std::vector<double> all_means(n_means * nprocs);

    int n = 0;
    if (process_no == 0) {
        std::ifstream input("numbers", std::ios::binary);
        n = file_size(input);
        std::vector<unsigned char> numbers_file(n);
        
        input.read(reinterpret_cast<char*>(&numbers_file[0]), n);
        input.close();
        if (nprocs > n) {
            std::cerr << "Too many processes" << std::endl;
            MPI_Finalize();
            return 1;
        } else if(nprocs < n){
            numbers_file.resize(nprocs);
        }

        //copy numbers_file to numbers
        for (int i = 0; i < n; i++) {
            numbers[i] = (int)numbers_file[i];
        }

        //print number_file
        std::cout << "Numbers: ";
        for (int i = 0; i < n; i++) {
            std::cout << numbers[i] << " ";
        }
        std::cout << std::endl;
        printf("Initialized centers: ");
        for (int i = 0; i < n_means; i++) {
            means[i] = numbers_file[i];
            printf("%d ", int(means[i]));
        }
        std::cout << std::endl;
    }

    // State: rank 0 has all numbers and means

    // Distibute one number to each process
    int number = 0;
    MPI_Scatter(&numbers[0], 1, MPI_INT, &number, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //std::cout << "Process " << process_no << " got number " << number << std::endl;
    
    // Broadcast means to all processes
    MPI_Bcast(&means[0], n_means, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // State: all processes have means and one number

    bool converged = false;
    while(!converged) {
        // Find nearest cluster
        int cluster_index = find_nearest_cluster(number, means);

        // recompute means ofall clusters

        int local_size[n_means] = {0};
        int local_sum[n_means] = {0};

        local_size[cluster_index] = 1;
        local_sum[cluster_index] = number;

        int cluster_sizes[n_means] = {0};
        int cluster_sums[n_means] = {0};

        MPI_Allreduce(&local_size[0], &cluster_sizes[0], n_means, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_sum[0], &cluster_sums[0], n_means, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        // print
        // std::cout << "Process " << process_no << " cluster_sizes: ";
        // for (int i = 0; i < n_means; i++) {
        //     std::cout << cluster_sizes[i] << " ";
        // }
        // std::cout << std::endl;

        // std::cout << "Process " << process_no << " cluster_sums: ";
        // for (int i = 0; i < n_means; i++) {
        //     std::cout << cluster_sums[i] << " ";
        // }
        // std::cout << std::endl;

        std::vector<double> new_means(n_means);

        for (int i = 0; i < n_means; i++) {
            new_means[i] = cluster_sums[i] / cluster_sizes[i];
        }

        // Check if converged
        converged = true;
        for (int i = 0; i < n_means; i++) {
            if (new_means[i] != means[i]) {
                converged = false;
                break;
            }
        }

        means = new_means;
    }

    //std::cout << "Process " << process_no << " converged" << std::endl;

    // rank 0 will pretty print the clusters
    if (process_no == 0) {
        // get clusters
        std::vector<std::vector<int>> clusters(n_means);
        for (int i = 0; i < n; i++) {
            int cluster_index = find_nearest_cluster(numbers[i], means);
            clusters[cluster_index].push_back(numbers[i]);
        }
        for (int i = 0; i < n_means; i++) {
            print_cluster(means[i], clusters[i]);
        }
    }

    MPI_Finalize();
    return 0;
}




