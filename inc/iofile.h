#ifndef READER_H
#define READER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

class iofile {
private:
    std::vector<std::string> data;  // Store file content
    bool file_opened;  // Status flag
public:
    iofile();  // Constructor
    bool readFile(const string& filename);  // Function to read file
    void output_data(const string& filename, const vector<vector<double>>& dt);  // Function to output data
    std::vector<std::string> getData() const;  // Function to get data
};

#endif
