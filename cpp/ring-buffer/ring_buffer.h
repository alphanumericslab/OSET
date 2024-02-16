#ifndef RINGBUFFER_H
#define RINGBUFFER_H

#include <vector>
#include <mutex>

class RingBuffer {
public:
    RingBuffer(size_t capacity, size_t inRate, size_t outRate, size_t inBlockSize, size_t outBlockSize, size_t outputOverlap);
    bool write(const std::vector<float>& data);
    bool read(std::vector<float>& data);

private:
    std::vector<float> buffer;
    size_t inIndex;
    size_t outIndex;
    size_t inRate;
    size_t outRate;
    size_t inBlockSize;
    size_t outBlockSize;
    size_t capacity;
    size_t outputOverlap;
    std::mutex mutex;
};

#endif // RINGBUFFER_H
