#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <cstdlib>

using namespace std;

#include "resampler.h"

// Function to print help message
void printHelp() {
    std::cout << "Usage: resampler_test [options]" << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -h, --help      Print this help message" << std::endl;
}

int main(int argc, char* argv[]) {
    // Check for command-line options
    if (argc > 1) {
        std::string option = argv[1];
        if (option == "-h" || option == "--help") {
            // Print the help message and exit
            printHelp();
            return 0;
        }
        else {
            std::cerr << "Error: Invalid option. Use -h or --help for usage information." << std::endl;
            return 1;
        }
    }

    // Sample test of the Resampler class
    Resampler resampler;
    const int inputLength = 100;
    double inputSignal[inputLength];
    double outputSignal[inputLength * 2]; // Resample with a ConversionRate of 2

    // Create a simple input signal (sine wave)
    for (int i = 0; i < inputLength; i++) {
        double t = double(i) / 100.0; // Time in seconds
        inputSignal[i] = sin(2 * M_PI * 2 * t); // 2 Hz sine wave
    }

    // Perform resampling
    double conversionRate = 2.0; // Resample to twice the input sample rate
    int order = 4; // Use fourth-order interpolation
    resampler.resample(inputSignal, outputSignal, inputLength, conversionRate, order);

    // Print the results
    std::cout << "Input Signal: ";
    for (int i = 0; i < inputLength; i++) {
        std::cout << inputSignal[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Resampled Signal: ";
    for (int i = 0; i < inputLength * 2; i++) {
        std::cout << outputSignal[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}
