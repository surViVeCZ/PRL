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
        double mean = 0;
        for (int j = 0; j < clusters[i].size(); j++) {
            mean += clusters[i][j];
        }
        mean /= clusters[i].size();
        means[i] = mean;
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


int main(int argc, char *argv[])
{
    int process_no;
    int nprocs;
    //create list of means of size n_means
    std::vector<double> means(n_means);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &process_no);
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
    int n = 0;
    if (process_no == 0) {
        n = file_size(input);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    int chunk_size = n / nprocs;
    if (process_no < n % nprocs) {
        chunk_size++;
    }

    // numbers of processes must the same as numbers in file
    if (nprocs > n) {
        std::cerr << "Too many processes" << std::endl;
        MPI_Finalize();
        return 1;
    }
    if (n > nprocs) {
        input.seekg(n-nprocs, ios::beg);
        n = nprocs;
    }
    //read input numbers from binary file
    std::vector<unsigned char> numbers_file(n);
    if (process_no == 0) {
        // Read in the data for the root process from the binary file
        input.read(reinterpret_cast<char*>(&numbers_file[0]), n);
        input.close();

         // Split numbers_file into 4 clusters
        std::vector<std::vector<unsigned char>> clusters(n_means);
        for (int i = 0; i < n; i++) {
            clusters[i % n_means].push_back(numbers_file[i]);
        }
     
        for (int i = 0; i < n_means; i++) {
            //print cluster mean from means list
            double mean = 0;
            for (int j = 0; j < clusters[i].size(); j++) {
                mean += clusters[i][j];
            }
            mean /= clusters[i].size();
            means[i] = mean;

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
    }
    // broadcast means to all processes
    MPI_Bcast(&means[0], n_means, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // assign each number to its nearest cluster
    std::vector<int> assignments(n);
    for (int i = 0; i < n; i++) {
        int cluster_index = find_nearest_cluster(numbers_file[i], means);
        assignments[i] = cluster_index;
    }

    if (process_no == 0) {
        std::cout << std::endl;
    }
    //call new_clusters until assignments dont change
    std::vector<int> new_assignments = new_clusters(process_no, n, numbers_file, assignments);
    while (assignments != new_assignments) {
        assignments = new_assignments;
        new_assignments = new_clusters(process_no, n, numbers_file, assignments);
    }
    new_assignments = new_clusters(process_no, n, numbers_file, assignments);
    MPI_Finalize();
    return 0;
}



