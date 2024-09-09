#include "ring_buffer.h"
#include <iostream>
#include <vector>

std::vector<float> generateRamp(size_t start, size_t size) {
    std::vector<float> ramp(size);
    for (size_t i = 0; i < size; ++i) {
        ramp[i] = static_cast<float>(start + i);
    }
    return ramp;
}

int main() {
    // Parameters for the ring buffer
    size_t capacity = 100;          // Total capacity of the ring buffer
    size_t inRate = 1;              // Input rate
    size_t outRate = 1;             // Output rate
    size_t inBlockSize = 10;        // Input block size
    size_t outBlockSize = 10;       // Output block size
    size_t outputOverlap = 2;       // Overlap size for output blocks

    // Number of blocks to write
    size_t numBlocks = 5;

    // Parameters for the ring buffer
    std::cout << "capacity = " << capacity << std::endl;
    std::cout << "inRate = " << inRate << std::endl;
    std::cout << "outRate = " << outRate << std::endl;
    std::cout << "inBlockSize = " << inBlockSize << std::endl;
    std::cout << "outBlockSize = " << outBlockSize << std::endl;
    std::cout << "outputOverlap = " << outputOverlap << std::endl;
 
   // Create the ring buffer
    RingBuffer ringBuffer(capacity, inRate, outRate, inBlockSize, outBlockSize, outputOverlap);

    // Write multiple blocks to the ring buffer
    std::cout << "Writing to ring buffer..." << std::endl;
    for (size_t i = 0; i < numBlocks; ++i) {
        std::vector<float> ramp = generateRamp(i * inBlockSize, inBlockSize);
        if (!ringBuffer.write(ramp)) {
            std::cerr << "Failed to write to ring buffer at block " << i << "." << std::endl;
            return -1;
        }
    }

    // Read and print data from the ring buffer
    std::cout << "Reading from ring buffer..." << std::endl;
    for (size_t i = 0; i < numBlocks; ++i) {
        std::vector<float> outputData;
        if (!ringBuffer.read(outputData)) {
            std::cerr << "Failed to read from ring buffer at block " << i << "." << std::endl;
            return -1;
        }

        // Print the output data
        std::cout << "Output Data Block " << i << ": ";
        for (size_t i = 0; i < outputData.size(); ++i) {
            std::cout << outputData[i] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}