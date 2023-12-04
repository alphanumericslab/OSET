// trimmed_filter_main.cpp for the trimmed_filter.cpp (for compiling and running from the terminal)

/*
% Revision History:
%   2008: First release
%   2023: Renamed from deprecated version TrimmedFilter.cpp
%
% Reza Sameni, 2022-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <cstring>

#include "list.h"
#include "trimmed_filter.h"

struct Arguments {
    std::string mode;
    int wlen;
    int alpha;
    std::string weightsFile;
    std::string inputFile;
    std::string outputFile;
};

bool parseArguments(int argc, char* argv[], Arguments& args) {
    if (argc < 2) {
        return false;
    }

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-mode") {
            if (i + 1 < argc) {
                args.mode = argv[i + 1];
                ++i;
            } else {
                return false;
            }
        } else if (arg == "-wlen") {
            if (i + 1 < argc) {
                args.wlen = std::atoi(argv[i + 1]);
                ++i;
            } else {
                return false;
            }
        } else if (arg == "-alpha") {
            if (i + 1 < argc) {
                args.alpha = std::atoi(argv[i + 1]);
                ++i;
            } else {
                return false;
            }
        } else if (arg == "-weights") {
            if (i + 1 < argc) {
                args.weightsFile = argv[i + 1];
                ++i;
            } else {
                return false;
            }
        } else if (arg == "-i") {
            if (i + 1 < argc) {
                args.inputFile = argv[i + 1];
                ++i;
            } else {
                return false;
            }
        } else if (arg == "-o") {
            if (i + 1 < argc) {
                args.outputFile = argv[i + 1];
                ++i;
            } else {
                return false;
            }
        }
    }

    return true;
}

int main(int argc, char* argv[]) {
    Arguments args;

    // Check if the user wants help
    if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0)) {
        std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
        std::cout << "Options:" << std::endl;
        std::cout << "  -h, --help     Display this help message" << std::endl;
        std::cout << "  -e, --example  Run the example with hardcoded input" << std::endl;
        std::cout << "  -mode <mode>   Filtering mode (mean, median, trmean, wmedian)" << std::endl;
        std::cout << "  -wlen <w>      Window length for filtering" << std::endl;
        std::cout << "  -alpha <alpha> Alpha parameter for trimmed mean (if applicable; used only in 'trmean' and 'wmedian' modes)" << std::endl;
        std::cout << "  -weights <weights_file> Weighting window file (if applicable)" << std::endl;
        std::cout << "  -i <input_file>         Input CSV file containing the time-series to be filtered" << std::endl;
        std::cout << "  -o <output_file>        Output CSV file containing the filtered time-series" << std::endl;
        std::cout << std::endl;
        std::cout << "Full Syntax:" << std::endl;
        std::cout << "  " << argv[0] << " -mode <mode> -wlen <w> [-alpha <alpha>] [-weights <weights_file>] -i <input_file> -o <output_file>" << std::endl;
        std::cout << std::endl;
        std::cout << "Examples:" << std::endl;
        std::cout << "  " << argv[0] << " -mode mean -wlen 5 -i input.csv -o output.csv" << std::endl;
        std::cout << "  " << argv[0] << " -mode median -wlen 3 -i input.csv -o output.csv" << std::endl;
        std::cout << "  " << argv[0] << " -mode trmean -wlen 5 -alpha 2 -i input.csv -o output.csv" << std::endl;
        std::cout << "  " << argv[0] << " -mode wmedian -wlen 5 -weights weights.csv -i input.csv -o output.csv" << std::endl;
        return 0;
    }

    // Check if the user wants to run the example
    if (argc == 2 && (strcmp(argv[1], "-e") == 0 || strcmp(argv[1], "--example") == 0)) {
        // Example data (hardcoded)
        double x[14] = {5.7, 7.5, 8, 8, 9, 17., -1, 22, 32, 1.5, 2, 9, -18, 7};
        int n = 14;
        double y[14];

        // Perform the filtering (example with trimmed mean)
        trimmed_filter(x, y, n, "trmean", 5, 2);

        // Print the filtered output
        std::cout << "Input:  Output:" << std::endl;
        for (int i = 0; i < n; ++i) {
            std::cout << x[i] << "\t" << y[i] << std::endl;
        }

        return 0;
    }

    // Parse command-line arguments
    if (!parseArguments(argc, argv, args)) {
        std::cout << "Invalid arguments. Use -h or --help for usage information." << std::endl;
        return 1;
    }

    // Check if the user provided all required options
    if (args.mode.empty() || args.inputFile.empty() || args.outputFile.empty()) {
        std::cout << "Error: Please provide all required options. Use '-h' or '--help' for usage details." << std::endl;
        return 1;
    }

    // Open the input CSV file and read the data
    std::ifstream infile(args.inputFile);
    if (!infile) {
        std::cerr << "Error: Failed to open input file: " << args.inputFile << std::endl;
        return 1;
    }

    // Read the input data from the CSV file into a vector
    std::vector<double> x;
    double value;
    while (infile >> value) {
        x.push_back(value);
    }
    infile.close();

    // Calculate the number of elements in the input vector
    int n = static_cast<int>(x.size());

    // Allocate memory for the output array
    double* y = new double[n];

    // Check if the filter mode is "wmedian" and if h values are provided
    double* h = nullptr; // Default value for h (nullptr if not used)

    if (args.mode == "wmedian") {
        if (!args.weightsFile.empty()) {
            // Open the weights CSV file and read the data
            std::ifstream hfile(args.weightsFile);
            if (!hfile) {
                std::cerr << "Error: Failed to open weights file: " << args.weightsFile << std::endl;
                delete[] y;
                return 1;
            }

            std::vector<double> h_values;
            while (hfile >> value) {
                h_values.push_back(value);
            }
            hfile.close();

            h = h_values.data();
        } else {
            if (args.wlen <= 0) {
                std::cerr << "Error: Weighting window size (w) must be greater than zero for 'wmedian' filter mode." << std::endl;
                delete[] y;
                return 1;
            }
        }
    }

    // Perform the filtering
    trimmed_filter(x.data(), y, n, args.mode.c_str(), args.wlen, args.alpha, h);

    // Open the output CSV file and write the filtered data
    std::ofstream outfile(args.outputFile);
    if (!outfile) {
        std::cerr << "Error: Failed to open output file: " << args.outputFile << std::endl;
        delete[] y;
        return 1;
    }

    // Write the filtered data to the output CSV file
    for (int i = 0; i < n; ++i) {
        outfile << y[i] << std::endl;
    }
    outfile.close();

    // Cleanup
    delete[] y;

    std::cout << "Filtering complete. Filtered data stored in " << args.outputFile << std::endl;

    return 0;
}
