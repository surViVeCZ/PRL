#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

//function calculates number of bytes of binary file
int file_size(ifstream& file){
    file.seekg(0, ios::end);
    int n = file.tellg();
    file.seekg(0, ios::beg);
    return n;
}

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    int process_no, N; //number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &process_no);
    MPI_Comm_size(MPI_COMM_WORLD, &N);

    //opening binary file
    std::ifstream input("numbers", std::ios::binary);
    
    // Check if the file was opened successfully
    if (!input){
        std::cerr << "Failed to open file" << std::endl;
        MPI_Finalize();
        return 1;
    }

    int n = file_size(input);
    int n_per_process = n / N;

  

    // Compute the median on the root process
    int median;
     //opening binary file
    if (process_no == 0) {
        std::ifstream input("numbers", std::ios::binary);
    
        // Check if the file was opened successfully
        if (!input){
            std::cerr << "Failed to open file" << std::endl;
            return 1;
        }

        int n = file_size(input);
        std::vector<unsigned char> numbers_file(n);
        input.read(reinterpret_cast<char*>(&numbers_file[0]), n);

        input.close();
        // Output the numbers to the console
        for (int i = 0; i < n; i++){
            std::cout << static_cast<int>(numbers_file[i]) << " ";
        }
        std::cout << std::endl;

        median = numbers_file[n/2];
        std::cout << "Median: " << median << std::endl;
    }
      // Determine the starting position and length of the data for this process
    int start_pos = process_no * n_per_process;
    int end_pos = (process_no == N - 1) ? n : (start_pos + n_per_process);
    int len = end_pos - start_pos;

    // Read in the data for this process from the binary file
    std::vector<unsigned char> numbers_file(len);
    input.seekg(start_pos);
    input.read(reinterpret_cast<char*>(&numbers_file[0]), len);
    input.close();

    //TODO rewrite using scatter
    //split numbers into same size groups and send 1 group to each process using MPI_Scatter (first 2 numbers to process 0, next 2 to process 1, etc.)
    // MPI_Scatter(&numbers_file[0], len, MPI_INT, &numbers_file[0], len, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast the median to all processes
    MPI_Bcast(&median, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Compute the L, E, and G groups for this process
    std::vector<unsigned char> L, E, G;
    for (int i = 0; i < len; i++) {
        if (numbers_file[i] < median) {
            L.push_back(numbers_file[i]);
        } else if (numbers_file[i] > median) {
            G.push_back(numbers_file[i]);
        } else {
            E.push_back(numbers_file[i]);
        }
    }

    // Print out the L, E, and G groups for this process
    std::cout << "Process " << process_no << ":" << std::endl;
    std::cout << "  L: ";
    for (auto num : L) {
        std::cout << static_cast<int>(num) << " ";
    }
    std::cout << std::endl;

    std::cout << "  E: ";
    for (auto num : E) {
        std::cout << static_cast<int>(num) << " ";
    }
    std::cout << std::endl;

    std::cout << "  G: ";
    for (auto num : G) {
        std::cout << static_cast<int>(num) << " ";
    }
    std::cout << std::endl;

    
    //compute sums of L, E, and G groups
    int sum_L = 0, sum_E = 0, sum_G = 0;
    for (auto num : L) {
        sum_L += num;
    }
    for (auto num : E) {
        sum_E += num;
    }
    for (auto num : G) {
        sum_G += num;
    }

    //merge sum of L of each process
    int sum_L_all = 0;
    int sum_E_all = 0;
    int sum_G_all = 0;
    MPI_Reduce(&sum_L, &sum_L_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_E, &sum_E_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_G, &sum_G_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);


    MPI_Finalize();
    return 0;
}
