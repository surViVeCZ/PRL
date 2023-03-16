#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>   
using namespace std;

//function calcultes number of bytes of binary file
int file_size(ifstream& file){
    file.seekg(0, ios::end);
    int length = file.tellg();
    file.seekg(0, ios::beg);
    return length;
}

int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    int process_no;
    int median;
    int num_processes;
    MPI_Comm_rank(MPI_COMM_WORLD, &process_no);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    int recvcounts[num_processes] = {1}; //number of elements to receive from each process - each sends only 1 elements
    int displs[num_processes]; //displacement of each element in the receive buffer
    for (int i = 0; i < num_processes; i++) {
        displs[i] = i;
    }

    //opening binary file
    if (process_no == 0) {
        ifstream input("numbers", ios::binary);
    
        // Check if the file was opened successfully
        if (!input){
            cerr << "Failed to open file" << endl;
            return 1;
        }

        
        int length = file_size(input);
        char buffer[64];
        input.read(buffer, 64);

        input.close();
        int numbers_file[length] = {0};
        // Output the numbers to the console
        for (int i = 0; i < length; i++){
            cout << static_cast<int>(buffer[i]) << " ";
            //store buffer[i] into int array
            numbers_file[i] = static_cast<int>(buffer[i]);
        }
        cout << endl;

        median = numbers_file[length/2];
        std::cout << "Median: " << median << std::endl;
    }

    // Broadcast the median to all processes
    MPI_Bcast(&median, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Perform some operation on the median (in this case, simply adding 1)
    median += 1;

    // Create an array to hold the updated median for each process
    int updated_medians[num_processes];
    // Gather the updated median from all processes to the root process
    MPI_Gatherv(&median, 1, MPI_INT, updated_medians, recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    if (process_no == 0) {
        std::cout << "Updated median: " << updated_medians[0] << std::endl;
    }

    MPI_Finalize();
    return 0;
}
