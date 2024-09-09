# Trimmed Filter

The `trimmed_filter` function is a C++ implementation of a filter that supports moving averaging, median filtering, trimmed mean filtering, and weighted median filtering. It is designed to process one-dimensional signals, providing a robust way to remove noise or outliers from data. See [Gonzalo R. Arce (2004)](https://www.amazon.com/Nonlinear-Signal-Processing-Statistical-Approach/dp/0471676241)
 for details.

## Function description

The `trimmed_filter` function is implemented in C++ and can be compiled using the provided `bash.sh` script or used as a MATLAB MEX function. It takes an input array representing the signal to be filtered and applies different nonlinear filtering modes based on the specified parameters.

### Supported filtering modes:

1. **Moving averaging (mean)**: Calculates the moving average of the signal using a specified sliding window length.

2. **Median filtering (median)**: Computes the median value within a specified sliding window.

3. **Trimmed mean filtering (trmean)**: Calculates the trimmed mean within a specified sliding window, excluding the first and last 'alpha' samples.

4. **Weighted median filtering (wmedian)**: Calculates the weighted median within a specified window using custom weights.

### Syntax:

```bash
./trimmed_filter -m [filtering_mode] -w [window_length] [-a alpha] [-i input_file] [-o output_file] [-weights weights_file]
```

### Arguments:

- `[filtering_mode]`: Specify the filtering mode. Use one of the following options: 'mean', 'median', 'trmean', 'wmedian'.

- `[window_length]`: The window length used for filtering (an integer value).

- `[alpha]`: The alpha parameter used for trimmed mean filtering (an integer value). Only applicable for 'trmean' mode.

- `[input_file]`: The CSV file containing the input signal. Specify this argument if you want to read the input from a file.

- `[output_file]`: The CSV file to store the filtered output.

- `[weights_file]`: The CSV file containing weights for weighted median filtering. Specify this argument only for 'wmedian' mode.

## Compiling with g++

To compile the C++ code using `g++`, navigate to the root folder of the project and run the following command in the terminal:

```bash
g++ trimmed_filter_main.cpp trimmed_filter.cpp list.cpp -o trimmed_filter
```

After successful compilation, you can use the `trimmed_filter` function from the terminal as shown in the example section below.

## Using MATLAB's MEX interface

To use the `trimmed_filter` function as a MATLAB MEX function, you need to compile it with the MATLAB MEX command. Open MATLAB and navigate to the root folder of the project. Then, run the following command in the MATLAB Command Window:

```matlab
mex trimmed_filter_mex.cpp trimmed_filter.cpp list.cpp -output trimmed_filter
```

After successful compilation, you can call the `trimmed_filter` function directly from MATLAB as shown in the example section below. Note that the `mex` function requires an external C/C++ compiler on the local system. It requires setup, if it has not been set up before:

```matlab
mex -setup
```


### Example usage from the terminal

Suppose we have a CSV file named `input_signal.csv` containing the input signal that we want to filter. The file looks like this:

```csv
5.7
7.5
8
8
9
17
-1
22
32
1.5
2
9
-18
7
```

To filter the input signal using trimmed mean mode, window length 5, and alpha 2, we can run the following command:

```bash
./trimmed_filter -m trmean -w 5 -a 2 -i input_signal.csv -o filtered_output.csv
```

The filtered output will be saved in a CSV file named `filtered_output.csv`.

Sure! Here's the MATLAB example that you can add to the README.md file:

### Example usage in MATLAB
This example filters a sample vector using the trimmed mean mode with window length 5 and alpha 2 (replace 'trimmed_filter' with the appropriate path and file name used to compile the code with `mex`):

```matlab
%% Example Usage in MATLAB

% Example data (same as in the README.md)
x = [5.7, 7.5, 8, 8, 9, 17, -1, 22, 32, 1.5, 2, 9, -18, 7];
wlen = 7;
alpha = 2;

% Call the trimmed_filter function in different modes
y_mean = trimmed_filter(x, 'mean', wlen);
y_median = trimmed_filter(x, 'median', wlen);
y_trmean = trimmed_filter(x, 'trmean', wlen, alpha);

window = hamming(wlen);
% remove the DC bias and normalize the window
window = window - window(1);
window = window/sum(window);
y_wmedian = trimmed_filter(x, 'wmedian', wlen, alpha, window);
y = [y_mean, y_median, y_trmean, y_wmedian];

% Display and plot the filtered output
disp('Input Signal:');
disp(x);
disp('Input and filtered outputs:');
disp([x', y]);

figure
plot(x);
hold on
plot(y);
grid
legend('input', 'mean', 'median', 'trmean', 'wmedian');
```

Remember to include the required steps to compile the MEX function as mentioned in the previous sections of the README.md file before running this MATLAB script.

## Revision history:

- 2008: First release
- 2023: Commented and renamed from deprecated version TrimmedFilter

## Further reading
Gonzalo R. Arce (2004). Nonlinear Signal Processing: A Statistical Approach. John Wiley & Sons Inc.


## Credit information:

Reza Sameni, 2008-2023 The Open-Source Electrophysiological Toolbox
https://github.com/alphanumericslab/OSET