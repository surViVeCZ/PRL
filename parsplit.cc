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

    // Determine the starting position and length of the data for this process
    int start_pos = process_no * n_per_process;
    int end_pos = (process_no == N - 1) ? n : (start_pos + n_per_process);
    int len = end_pos - start_pos;

    // Allocate memory for local data
    std::vector<unsigned char> local_numbers(len);

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
        MPI_Scatter(&numbers_file[0], len, MPI_UNSIGNED_CHAR,
                    &local_numbers[0], len, MPI_UNSIGNED_CHAR,
                    0, MPI_COMM_WORLD);
    } else {
        // Allocate memory for root process data
        std::vector<unsigned char> numbers_file(n);

        // Scatter the data to all processes
        MPI_Scatter(NULL, len, MPI_UNSIGNED_CHAR,
                    &local_numbers[0], len, MPI_UNSIGNED_CHAR,
                    0, MPI_COMM_WORLD);
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
    for (int i = 0; i < len; i++) {
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
    

    // Gather the L, E, and G groups from all processes
    std::vector<unsigned char> L_all, E_all, G_all;
    L_all.resize(N);
    E_all.resize(N);
    G_all.resize(N);

    // MPI_Gather(&L[0], L.size(), MPI_UNSIGNED_CHAR,
    //            &L_all[0], L.size(), MPI_UNSIGNED_CHAR,
    //            0, MPI_COMM_WORLD);
    
    // MPI_Gather(&E[0], E.size(), MPI_UNSIGNED_CHAR,
    //             &E_all[0], E.size(), MPI_UNSIGNED_CHAR,
    //             0, MPI_COMM_WORLD);
    
    // MPI_Gather(&G[0], G.size(), MPI_UNSIGNED_CHAR,
    //             &G_all[0], G.size(), MPI_UNSIGNED_CHAR,
    //             0, MPI_COMM_WORLD);


    // // Merge the L, E, and G groups into L_all, E_all, and G_all on the root process
    // if (process_no == 0) {
    //     // Determine the total size of L, E, and G across all processes
    //     int total_L_size, total_E_size, total_G_size;
    //     MPI_Reduce(&L.size(), &total_L_size, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    //     MPI_Reduce(&E.size(), &total_E_size, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    //     MPI_Reduce(&G.size(), &total_G_size, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    //     // Resize L_all, E_all, and G_all to the total size of L, E, and G
    //     L_all.resize(total_L_size);
    //     E_all.resize(total_E_size);
    //     G_all.resize(total_G_size);

    //     // Gather the L, E, and G groups from all processes again, this time into L_all, E_all, and G_all
    //     MPI_Gather(&L[0], L.size(), MPI_UNSIGNED_CHAR,
    //             &L_all[0], L.size(), MPI_UNSIGNED_CHAR,
    //             0, MPI_COMM_WORLD);
    //     MPI_Gather(&E[0], E.size(), MPI_UNSIGNED_CHAR,
    //             &E_all[0], E.size(), MPI_UNSIGNED_CHAR,
    //             0, MPI_COMM_WORLD);
    //     MPI_Gather(&G[0], G.size(), MPI_UNSIGNED_CHAR,
    //             &G_all[0], G.size(), MPI_UNSIGNED_CHAR,
    //             0, MPI_COMM_WORLD);
    // }




    MPI_Finalize();
    return 0;
}
