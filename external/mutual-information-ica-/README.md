#  MILCA and SNICA: Independent Component Analysis Algorithms

**COPYRIGHT NOTICE**: *This directory includes files, programs, and data that are sourced from the public domain for ACADEMIC PURPOSES (see references and list of contributors below). **These files are excluded from the OSET open-source license and fall under their respective ownership and permissions**. Users are advised to contact the original owners of these codes on a case-by-case basis to obtain permission for reproducing or utilizing them. It is also important to appropriately acknowledge the original developers in any usage or publications that rely on the files in this directory.*

**NOTE**: *The current codebase has been downloaded in 2006 from [here](https://www.ucl.ac.uk/ion/milca-0). The source codes might be out-of-date; consider using more recent versions.*

## Overview
MILCA (Mutual Information Least-dependent Component Analysis) and SNICA (Stochastic Non-negative Independent Component Analysis) are Independent Component Analysis (ICA) algorithms. They leverage an accurate Mutual Information (MI) estimator to identify the least dependent components under a linear transformation. SNICA uses a non-negativity constraint.

The MI estimator is data-efficient, adaptive, and has minimal bias. This estimator enables contrast functions for ICA and can be used in performance and reliability tests, as well as for cluster analysis. The package provides the algorithms for these operations and extends to include time structure, often found in physical data. For more details, refer to the references and the contributors' publications.

## Codebase

- **[MILCA (standard ICA)](./milca.C)**: uses only the instantaneous information in the signal, optimal for signals with white spectra (no time structure)
  - **Input**: Original (mixed) components
  - **Output**: Least dependent components, de-mixing matrix
  
- **[MILCAdelay](./milcadelay.C)**: uses in addition any time structure in the signals. Can separate also two Gaussians with different spectra [1]
  - **Input**: Original (mixed) components
  - **Output**: Least dependent components, de-mixing matrix
  
- **[ICATests (Reliability Tests)](./ICAtests.C)**: can be used for any ICA-output, not necessarily from MILCA [2],[1]
  - **Input**: Any ICA-output
  - **Output**: Dependency matrix, variability matrix, MI vs rotation angle plots
  
- **[MIClustering](./MIClustering.C)**: hierarchical clustering algorithm based on the grouping property of MI [4],[1]
  - **Input**: Any ICA-output
  - **Output**: Dendrogram of dependencies

- **[MIxnyn](./MIxnyn.C)**: calculates mutual information between two input channels of arbitrary dimensions [3],[1]
  - **Input**: Multidimensional input signal
  - **Output**: MI value

- **[MIhigherdim](./MIhigherdim.C)**: calculates mutual information (redundancy) between any number of one-dimensional input channels [3],[1]
  - **Input**: Multidimensional input signal
  - **Output**: MI value

## Installation

To compile the code, you'll need a C++ compiler. The project includes a Makefile to simplify the build process.

### macOS/Linux

1. Open a terminal and navigate to the directory containing the Makefile.
2. Run the following command to compile the package:
    ```bash
    make
    ```

### Windows

For Windows, you may use MinGW or any other C++ compiler compatible with GNU Make.

1. Open Command Prompt and navigate to the directory containing the Makefile.
2. Run the following command to compile the package:
    ```bash
    mingw32-make
    ```

If you don't have GNU Make installed, you can manually compile each source file using `g++` as shown in the Makefile.

### Executables

After successfully running the `make` or `mingw32-make` command, the following executables will be generated:

- `milca`
- `milcadelay`
- `ICAtests`
- `MIhigherdim`
- `MIxnyn`
- `MIClustering`


For more detailed instructions on generating the executables, you can refer to the [original codebase](https://www.ucl.ac.uk/ion/milca-0).

## References
1- Least Dependent Component Analysis Based on Mutual Information Harald Stögbauer, Alexander Kraskov, Sergey A. Astakhov, and Peter Grassberger , Phys. Rev. E 70 (6)  066123, 2004

2- Reliability of ICA estimates with mutual information  Harald Stögbauer, Ralph G. Andrzejak, Alexander Kraskov and Peter Grassberger, Independent Component Analysis and Blind Signal Separation Lecture Notes in Computer Science 3195: 209-216, 2004

3- Estimating mutual information A. Kraskov, H. Stögbauer, and P. Grassberger,  Phys. Rev. E 69 (6) 066138, 2004

4- Hierarchical clustering using mutual information A. Kraskov, H. Stögbauer, R. G. Andrzejak, and P. Grassberger , Europhysics Letters 70 (2): 278-284,  2005

5- Spectral Mixture Decomposition by Least Dependent Component Analysis Sergey A. Astakhov, Harald Stögbauer, Alexander Kraskov,  and Peter Grassberger

6- Monte Carlo Algorithm for Least Dependent Non-Negative Mixture Decomposition Sergey A. Astakhov, Alexander Kraskov, Harald Stögbauer, and P. Grassberger , Analytical Chemistry Web Release Date: 19-Jan-2006

## FOR SCIENTIFIC USE ONLY
These codes are free for research and education purposes. Commercial or military use is prohibited.

## NO WARRANTY
The software is provided "as-is," without any express or implied warranty. We are not liable for any damages arising from its use.

---

## Contributors

- Sergey Astakhov
- Peter Grassberger
- Alexander Kraskov
- Harald Stögbauer