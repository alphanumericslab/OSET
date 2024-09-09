#include "ring_buffer.h"
#include <iostream>

RingBuffer::RingBuffer(size_t capacity, size_t inRate, size_t outRate, size_t inBlockSize, size_t outBlockSize, size_t outputOverlap)
    : capacity(capacity), inRate(inRate), outRate(outRate), inBlockSize(inBlockSize), 
        outBlockSize(outBlockSize), outputOverlap(outputOverlap), inIndex(0), outIndex(0) {
    buffer.resize(capacity, 0.0f);
}


bool RingBuffer::write(const std::vector<float>& data) {
    std::lock_guard<std::mutex> lock(mutex);
    
    if (data.size() != inBlockSize) {
        std::cerr << "Input block size mismatch." << std::endl;
        return false;
    }

    for (size_t i = 0; i < inBlockSize; ++i) {
        buffer[inIndex] = data[i];
        inIndex = (inIndex + 1) % capacity;
        
        if (inIndex == outIndex) {
            std::cerr << "Buffer overflow." << std::endl;
            return false;
        }
    }

    return true;
}

bool RingBuffer::read(std::vector<float>& data) {
    std::lock_guard<std::mutex> lock(mutex);
    data.clear();

    size_t availableSamples = (inIndex + capacity - outIndex) % capacity;
    size_t expectedSamples = outBlockSize + outputOverlap;

    if (availableSamples < expectedSamples) {
        std::cerr << "Not enough data available." << std::endl;
        return false;
    }

    for (size_t i = 0; i < expectedSamples; ++i) {
        data.push_back(buffer[outIndex]);
        outIndex = (outIndex + 1) % capacity;
    }

    // Adjust for output rate
    outIndex = (outIndex + capacity - (outBlockSize * outRate / inRate)) % capacity;
    
    return true;
}

