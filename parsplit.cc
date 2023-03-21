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

    // Allocate memory for local data
    std::vector<unsigned char> local_numbers(n_per_process);

    if (process_no == 0) {
        // Read in the data for the root process from the binary file
        std::vector<unsigned char> numbers_file(n);
        input.read(reinterpret_cast<char*>(&numbers_file[0]), n);
        input.close();

        // Output the numbers to the console
        for (int i = 0; i < n; i++){
            std::cout << static_cast<int>(numbers_file[i]) << " ";
        }
        std::cout << std::endl;

        // Scatter the data to all processes
        MPI_Scatter(&numbers_file[0], n_per_process, MPI_UNSIGNED_CHAR,
                    &local_numbers[0], n_per_process, MPI_UNSIGNED_CHAR,
                    0, MPI_COMM_WORLD);
    } else {
        MPI_Scatter(0, n_per_process, MPI_UNSIGNED_CHAR,&local_numbers[0], n_per_process, MPI_UNSIGNED_CHAR,0, MPI_COMM_WORLD);
    }

    // Compute the median on the root process
    int median;
    if (process_no == 0) {
        median = local_numbers[n_per_process/2];
        std::cout << "Median: " << median << std::endl;
    }

    // Broadcast the median to all processes
    MPI_Bcast(&median, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Compute the L, E, and G groups for this process
    std::vector<unsigned char> L, E, G;
    for (int i = 0; i < n_per_process; i++) {
        //print local_numbers
        if (local_numbers[i] < median) {
            L.push_back(local_numbers[i]);
        } else if (local_numbers[i] == median) {
            E.push_back(local_numbers[i]);
        } else {
            G.push_back(local_numbers[i]);
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
    

    // Allocate a new buffer to hold all the gathered data on the root process
    std::vector<unsigned char> all_numbers(n);
    std::vector<unsigned char> L_all, E_all, G_all;

    // Gather the data from all the processes into the all_numbers buffer on the root process
    MPI_Gather(&local_numbers[0], n_per_process, MPI_UNSIGNED_CHAR,
            &all_numbers[0], n_per_process, MPI_UNSIGNED_CHAR,
            0, MPI_COMM_WORLD);
    // Merge the E, G, and L groups on the root process
    if (process_no == 0) {
        for (int i = 0; i < n; i++) {
            if (all_numbers[i] < median) {
                L_all.push_back(all_numbers[i]);
            } else if (all_numbers[i] == median) {
                E_all.push_back(all_numbers[i]);
            } else {
                G_all.push_back(all_numbers[i]);
            }
        }
    }

    // Print out the E_all, G_all, and L_all groups in the root process
    if (process_no == 0) {
        std::cout << std::endl;
        std::cout << "L_all: ";
        for (auto num : L_all) {
            std::cout << static_cast<int>(num) << " ";
        }
        std::cout << std::endl;

        std::cout << "E_all: ";
        for (auto num : E_all) {
            std::cout << static_cast<int>(num) << " ";
        }
        std::cout << std::endl;

        std::cout << "G_all: ";
        for (auto num : G_all) {
            std::cout << static_cast<int>(num) << " ";
        }
        std::cout << std::endl;
    }

    MPI_Finalize();
    return 0;
}
