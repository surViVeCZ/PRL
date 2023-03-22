#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;


//function to merge L_all, E_all, G_all into root process and print it
void merge_and_print(vector<unsigned char> &L_all, vector<unsigned char> &E_all, vector<unsigned char> &G_all, int N, int n_per_process, int n){
    //resize buffer
    vector<unsigned char> all_numbers(n);
    //merge L, E, G groups
    merge(L_all.begin(), L_all.end(), E_all.begin(), E_all.end(), all_numbers.begin());
    merge(all_numbers.begin(), all_numbers.end(), G_all.begin(), G_all.end(), all_numbers.begin());
    //print merged vector
    cout << "Merged vector: ";
    for(int i = 0; i < n; i++){
        cout << static_cast<int>(all_numbers[i]) << " ";
    }
    cout << endl;
}

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
        MPI_Scatter(nullptr, n_per_process, MPI_UNSIGNED_CHAR,&local_numbers[0], n_per_process, MPI_UNSIGNED_CHAR,0, MPI_COMM_WORLD);
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
    std::vector<unsigned char> L_all;
    std::vector<unsigned char> E_all;
    std::vector<unsigned char> G_all;


    int L_size = L.size();
    int E_size = E.size();
    int G_size = G.size();
    std::vector<int> L_sizes(N);
    std::vector<int> E_sizes(N);
    std::vector<int> G_sizes(N);

    //gathering sizes of L, E, G groups
    MPI_Gather(&L_size, 1, MPI_INT, &L_sizes[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&E_size, 1, MPI_INT, &E_sizes[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&G_size, 1, MPI_INT, &G_sizes[0], 1, MPI_INT, 0, MPI_COMM_WORLD);

    //displacement of L, E, G groups in the buffer
    std::vector<int> L_offset(N);
    std::vector<int> E_offset(N);
    std::vector<int> G_offset(N);
    //specifies starting position of each group in the buffer
    for(int i = 1; i < N; i++){
        L_offset[i] = L_offset[i-1] + L_sizes[i-1];
        E_offset[i] = E_offset[i-1] + E_sizes[i-1];
        G_offset[i] = G_offset[i-1] + G_sizes[i-1];
    }

    //now we gather L, E, G groups using gatherv
    if(process_no == 0){
        //resize buffers
        L_all.resize(L_offset[N-1] + L_sizes[N-1]);
        E_all.resize(E_offset[N-1] + E_sizes[N-1]);
        G_all.resize(G_offset[N-1] + G_sizes[N-1]);

        MPI_Gatherv(&L[0], L_size, MPI_UNSIGNED_CHAR, &L_all[0], &L_sizes[0], &L_offset[0], MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
        MPI_Gatherv(&E[0], E_size, MPI_UNSIGNED_CHAR, &E_all[0], &E_sizes[0], &E_offset[0], MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
        MPI_Gatherv(&G[0], G_size, MPI_UNSIGNED_CHAR, &G_all[0], &G_sizes[0], &G_offset[0], MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    }
    else{
        L_all.resize(0);
        E_all.resize(0);
        G_all.resize(0);
        MPI_Gatherv(&L[0], L_size, MPI_UNSIGNED_CHAR, nullptr, nullptr, nullptr, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
        MPI_Gatherv(&E[0], E_size, MPI_UNSIGNED_CHAR, nullptr, nullptr, nullptr, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
        MPI_Gatherv(&G[0], G_size, MPI_UNSIGNED_CHAR, nullptr, nullptr, nullptr, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
    }

    // Print out the E_all, G_all, and L_all groups in the root process
    if (process_no == 0) {
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

        std::cout << "L_all: ";
        for (auto num : L_all) {
            std::cout << static_cast<int>(num) << " ";
        }
        std::cout << std::endl;
    }

    //initialize numbers_solution vector
    std::vector<unsigned char> numbers_solution(n);
    //get size L_all + E_all
    int L_all_size = L_all.size();
    int E_all_size = E_all.size();
  
    //merge all 3 groups (L_all, E_all, G_all) to numbers_solution
    std::merge(L_all.begin(), L_all.end(), E_all.begin(), E_all.end(), numbers_solution.begin());
    std::merge(numbers_solution.begin(), numbers_solution.begin() + L_all_size + E_all_size, G_all.begin(), G_all.end(), numbers_solution.begin());

    //print numbers_solution
    if(process_no == 0){
        std::cout << std::endl;
        std::cout << "numbers_solution: ";
        for (auto num : numbers_solution) {
            std::cout << static_cast<int>(num) << " ";
        }
        std::cout << std::endl;
    }

    MPI_Finalize();
    return 0;
}
