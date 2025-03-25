#include "iofile.h"

iofile::iofile() { }

bool iofile::readFile(const string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        data.push_back(line);
    }
    file.close();
    file_opened = true;
    return true;
}

std::vector<std::string> iofile::getData() const
{
    return data;
}

void iofile::output_data(const string& filename, const vector<vector<double>>& dt)
{
    std::ofstream file(filename);  // Open file for writing (overwrite mode)
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
    } else {
        for (auto& ee: dt)
        {
            for (auto& uu: ee)
            {
                file << to_string(uu) << " ";
            }
            file << std::endl;

        }
    }

    file.close();
}