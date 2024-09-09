# Fractional Resampler

The Resampler is a C++ class that provides functionality for resampling continuous signals with different interpolation orders. It can be used to convert signals from one sample rate to another using various interpolation techniques, such as Second Order (Triangular Convolution), Third Order (Quadratic Convolution), and Fourth Order. The function supports fractional sampling (non-integer ratios).

## Introduction

The Resampler class allows for efficient resampling of continuous signals. It can be used to process streaming data block by block, making it suitable for real-time applications or working with large datasets, where using standard resampling methods (upsampling/downsampling, followed by anti-imaging/anti-aliasing filtering) would be inefficient due to memory limitations.

The Resampler class has been tested and verified with C++ compilers such as g++ and clang++.

## Features

- Supports different interpolation orders (2nd to 4th order).
- Efficiently resamples continuous signals block by block.
- Simple API for easy integration into other applications.
- Supports both double precision and integer input signals.
- Provides a sample test program to demonstrate the resampling functionality.

## Files

The Resampler package contains the following files:

1. `resampler.h`: The header file containing the Resampler class declaration.
2. `resampler.cpp`: The implementation file for the Resampler class.
3. `resampler_test.cpp`: A sample test program demonstrating how to use the Resampler class.
4. `resampler_mex.cpp`: The MATLAB mexFunction interface for using the Resampler class in MATLAB.

## Usage

To use the Resampler class in your own project, include the `resampler.h` header file in your source code. Create an instance of the Resampler class and call the `resample` function to perform the resampling.

The `resample` function has the following signature:

```cpp
int resample(double* x, double* y, int InputLength, double ConversionRate, int Order);
```

- `x`: Pointer to the input signal array (double precision).
- `y`: Pointer to the output signal array (resampled signal).
- `InputLength`: Length of the input signal `x`.
- `ConversionRate`: Conversion rate for resampling (input sampling rate divided by the output sampling rate), i.e., ConversionRate is the fraction Ts(out)/Ts(in) or Fs(in)/Fs(out).
- `Order`: Order of interpolation (2: Second Order, 3: Third Order, 4: Fourth Order, default is 2).

## Compilation

To compile the sample test program (`resampler_test.cpp`) and the MATLAB mexFunction interface (`resampler_mex.cpp`) along with the Resampler class, use a C++ compiler such as g++ or clang++. Follow the steps below:

### Compiling the Sample Test Program (resampler_test.cpp):

1. Open a terminal or command prompt.
2. Navigate to the directory containing the Resampler files (`resampler.h`, `resampler.cpp`, and `resampler_test.cpp`).
3. Use the following command to compile the program:

   ```bash
   g++ resampler_test.cpp resampler.cpp -o resampler_test
   ```

   or

   ```bash
   clang++ resampler_test.cpp resampler.cpp -o resampler_test
   ```

4. After successful compilation, you will find an executable named `resampler_test` in the same directory.

### Compiling the MATLAB mexFunction Interface (resampler_mex.cpp):

1. Open MATLAB.
2. Navigate to the directory containing the Resampler files (`resampler.h`, `resampler.cpp`, and `resampler_mex.cpp`).
3. Use the following command to compile the `mexFunction` interface:

   ```matlab
   mex resampler_mex.cpp resampler.cpp
   ```

   This command will compile the C++ code into a MEX-file (`resampler_mex.mexa64` on Linux, `resampler_mex.mexmaci64` on macOS, or `resampler_mex.mexw64` on Windows).
   
   The compiled library name can be cnaged to an arbitrary name using the `-output` option in `mex`, e.g, 
      ```matlab
   mex resampler_mex.cpp resampler.cpp -output fractional_resampler
   ```
   which will change the compiled library names to `fractional_resampler.mexa64` on Linux, `fractional_resampler.mexmaci64` on macOS, or `fractional_resampler.mexw64` on Windows.

## MATLAB Usage

After compiling the `resampler_mex.cpp` with `mex`, you can use the Resampler class in MATLAB:
   ```matlab
   % Sample test of the Resampler class using the mexFunction interface
   inputSignal = sin(2 * pi * 2 * (0:0.01:1)); % 2 Hz sine wave
   conversionRate = 2.0; % Resample to half the input sample rate (input/output sampling frequency)
   order = 4; % Use fourth-order interpolation

   % Call the resampler_mex function to perform resampling
   outputSignal = resampler_mex(inputSignal, conversionRate, order);

   % Plot the original and resampled signals
   tInput = linspace(0, 1, length(inputSignal));
   tOutput = linspace(0, 1, length(outputSignal));
   plot(tInput, inputSignal, 'b-', tOutput, outputSignal, 'r--');
   legend('Input Signal', 'Resampled Signal');
   ```

Please note that the `resampler_test` program provides a more comprehensive test of the Resampler class and can be executed directly from the terminal as described in the [Compilation](#compilation) section.

---

## Revision history:

- 2012: First release
- 2023: Commented and renamed from deprecated version resampler.cpp

## Credit information:

Reza Sameni, 2012-2023 The Open-Source Electrophysiological Toolbox
https://github.com/alphanumericslab/OSET

---