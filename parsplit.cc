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
    MPI_Comm_rank(MPI_COMM_WORLD, &process_no);
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

        std::cout << "Median: " << numbers_file[length/2] << std::endl;
    }

    std::cout << "Hello from process " << process_no << "!" << std::endl;

    MPI_Finalize();
    return 0;
}
