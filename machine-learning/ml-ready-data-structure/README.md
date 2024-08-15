# A Generic Data Structure for AI-ML-Ready Feature Tables

## Overview
Software version control is crucial in both software engineering and data science. However, traditional version control methods often fall short in linking data and code and tracking the exact parameters used for training models and extracting features in artificial intelligence (AI) and machine learning (ML) applications. In many applications, even minor alterations in model parameters can lead to different and irreproducible outcomes. While source code can be effectively managed and tracked, the precise parameters and datasets used for model evaluation can change, often subtly but significantly, impacting results.

To ensure the reproducibility of feature extraction and model training, it is crucial to combine the exact software version, the exact parameters used to run the code, and the dataset on which the software is applied.

We propose a lightweight and generic data structure for storing AI-ML-ready flat feature tables. This structure is primarily designed for extracting features and generating flat feature tables from physiological time-series data but is applicable to a variety of other use cases. Our proposed data structure integrates codebase hashs, data record hashes, the parameters used to run the codes, and the interdependencies between methods executed sequentially on the record to produce the output feature tables.

## Data structure
```
analysis
├── record[1]
│   ├── name: record name (string)
│   ├── rel_path: relative path (string)
│   ├── md5chsum: record file MD5 checksum (HEX string)
│   ├── sampling_freq: record sampling frequency (float)
│   ├── mains_freq: powerline frequency (float)
│   ├── num_ch: number of channels (int)
│   └── methods: a list of methods/functions
│       ├── method[1]
│       │   ├── name: method/function name (string)
│       │   ├── rel_path: method/function relative path to codebase root (string)
│       │   ├── codebase_md5chsum: codebase MD5 checksum (HEX string)
│       │   ├── codebase_git_repo: Git repository URL/link (string)
│       │   ├── codebase_git_commit_id: Unique git commit ID (HEX string)
│       │   ├── inputs: list of inputs
│       │   │   ├── input[1]: 1st input
│       │   │   │   ├── name: input name (string). Can be one of these: 1) 'record' (refers to current record); 2) 'method_p_output_q' (output of method[p] output[q]), e.g., method_3_output_4
│       │   │   │   ├── type: input type (string): 'bool', 'int', 'float' or 'string'
│       │   │   │   └── channels: list of channels from the input (int array)
│       │   │   ├── ...
│       │   │   └── input[m]: m-th input
│       │   │       └── ...
│       │   ├── params
│       │   │   ├── param[1]: 1st parameter
│       │   │   │   ├── name: parameter name (string). Use internal name of the parameter
│       │   │   │   ├── type: parameter type (string): 'bool', 'int', 'float' or 'string'
│       │   │   │   └── value: parameter value (bool/int/float/string, single-valued or array)
│       │   │   ├── ...
│       │   │   ├── param[l]: l-th parameter
│       │   │   └── ...
│       │   ├── outputs: list of outputs
│       │   │   ├── output[1]: first output
│       │   │   │   ├── name: first output name. Use internal variable name of the output
│       │   │   │   ├── type: output type (string): 'bool', 'int', 'float' or 'string'
│       │   │   │   └── value: output value (bool/int/float/string, single-valued or array)
│       │   │   ├── ...
│       │   │   └── output[n]
│       │   ├── errors: collects warning/error messages returned by the method returned through standard IO(string or string array)
│       │   └── success: whether the method succeeded or failed (bool)
│       ├── method[2]
│       │   └── ...
│       ├── ...
│       └── method[k]
│           └── ...
├── ...
└── record[r]
    └── ...
```

## Implementations
Implementations of the AI-ML-ready data structure are available in [MATLAB](./matlab) and [Python](./python). 
After running the feature extraction codes in MATLAB and Python on the same database, you can merge the outputs to create a single JSON file that includes all the extracted features. Please run this [script](./python/merge_json_files.ipynb).

## Contributors
- [Somayyeh Mousavi](seyedeh.somayyeh.mousavi@emory.edu)
- [Reza Sameni](rsameni@dbmi.emory.edu)