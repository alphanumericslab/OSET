# Python-Based Feature Extraction for AI and ML Learning Applications

## Contents
This folder contains Python scripts for applying various feature extraction functions to a dataset. Additionally, using a simple example, we will explain how to modify the scripts based on your specific needs.

## Run the Codes
1. Clone the [OSET repository](https://github.com/alphanumericslab/OSET).
2. Run [main_python](main_python.py). You will see the JSON file in the output folder.
3. In the following table, you see the three main Python scripts of the ready feature extraction pipeline that you should modify to run on your dataset, based on the feature extraction functions that you have considered.
   
| File Name                   | Description                                                                                                    |
|-----------------------------|----------------------------------------------------------------------------------------------------------------|
| [main_python](main_python.py)              | The script initializes parameters, searches for data files, validates inputs, and calls a feature extraction function to process the data, then saves the results. |
|[extract_feature](extract_feature.py) | The script loads the data and then runs the feature extraction algorithms on it.|
|[function_callers](function_callers.py) |The script includes functions that call a specific function.| 

## Modify the Codes
1. Open the [main_python](python/main_python.py).
2. You should modify `Part 1: Setup and Initialization`. According to the aim of your project, you should define the following 4 parts:
   - Data_Input_Path: Path containing the data
   - Data_Output_Path: Path where the results will be saved
   - Desired_Output_File: Name of the output file
   - The given extension (for example, .mat, .csv, or .dat) to list the files in the folder and subfolders of Data_Input_Folder
3. Then, you should modify `Part 2: Parameter Configuration` to correct the sampling frequency, number of channels, and powerline frequency. Then, define the list of functions needed for feature extraction. `load` is a necessary function (Don't change it). Please add the other functions in the format provided below to the `Params['function_dict']`. All function inputs should be strings.

```
# Change n, the number that you assign to a function (e.g., 1, 2, 3, ...)
Fun_n: {
    'name': 'Fill this part, method/function name (string)', ...
    'rel_path': 'Fill this part, method/function relative path to codebase root (string)', ...
    'codebase_git_repo': 'Fill this part, Git repository URL/link (string)', ...
    'codebase_git_commit_id': 'Fill this part, Unique git commit ID (HEX string)' ...
}
```
**Please ensure that the function names are accurate. In the next stage they are called exactly as specified.**

`Fun_1` and `Fun_2` are functions related to example codes. You can deactivate or change them.
         
1. Open the [extract_feature](extract_feature.py).
1. You need to write scripts to call the functions and then save the outputs. Each function has four parts:
- Part 1: Function Setup and Initialization
```python
  # Part 1: Function Setup and Initialization 
  # Change n, the number related to the function that you have defined in main_python file . 
  method = list_functions['Fun_n']
  # Don't change 
  Num_fun += 1 
  feature['methods'][f'method_{Num_fun}'] = {
     'name': method['name'],
     'rel_path': method['rel_path'],
     'codebase_md5chsum': method['codebase_md5chsum'],
     'codebase_git_repo': method['codebase_git_repo'],
     'codebase_git_commit_id': method['codebase_git_commit_id']
     }
````
- Part 2: Call the function, choose the appropriate option based on whether you want to use whole data or process channels one by one (Whole_data, Channel_data)

```python
# Part 2: Whole_data
# Call the function to process the data
# Change n, the number related to the function that you have defined in main_python file.  
result_n = functioncaller_n(method['name'], data, file_name, params)

# Store the input info of the function
feature['methods'][f'method_{Num_fun}']['inputs'] = {
   # Don't change
   'name': 'whole_data',
   # Change
   'type': [Fill this part, input type (string): 'bool', 'int', 'float' or 'string'],
   # Don't change. We considered 0 for the list of all channels from the input
   'channels': 0
        }
```

```python
# Part 2: Channel_data
# Don't change 
channel_name = [None] * params['Number_Channels']
warning = [None] * params['Number_Channels']
success = [0] * params['Number_Channels']

# Add the line of code below, based on the number of outputs, and replace the placeholders with the internal names of the outputs
# ** output = [None] * params['Number_Channels'] **
output1 = [None] * params['Number_Channels']
output2 = [None] * params['Number_Channels']

# Process each channel individually
for k in range(params['Number_Channels']):
            
   # Call the function
   # Change n, the number related to the function that you have defined in main_python file.  
   result_n = functioncaller_n(method['name'], data, k, file_name, params)
   # Don't change
   channel_name[k] = k + 1
   # Change n, the number related to the function that you have defined in main_python file.  
   success[k] = result_n['success']
   warning[k] = result_n['warning']
   # Add the line of code below, based on the number of the outputs, and replace the placeholders with the internal names of the outputs
   # ** output1[k] = result_n['output1'] **
   # Change n
   output1[k] = result_n['output1']
   output2[k] = result_n['output2']
        
# Store the input info of the function
# Don't change
feature['methods'][f'method_{Num_fun}']['inputs'] = {
   # Don't change
   'name': 'channel_data',
   # Change
   'type': [Fill this part, input type (string): 'bool', 'int', 'float' or 'string'],
   # Don't change
   'channels': channel_name
        }
```

- Part 3: Save the parameters.

If there is one parameter, use the below structure.
```
# Store the params of the function
feature['methods'][f'method_{Num_fun}']['params'] = {
   'name': '[]', # Use internal name of the parameter,
   'type': '[]', # 'bool', 'int', 'float' or 'string'
   'val': []  # bool/int/float/string, single-valued or array
        }
```
If there are multiple parameters, use the below structure.
```
# Store the params of the function
feature['methods'][f'method_{Num_fun}']['params'] = {
   'param_1': {
      'name': '[]', # Use internal name of the parameter,
      'type': '[]', # 'bool', 'int', 'float' or 'string'
      'val': []  # bool/int/float/string, single-valued or array
            },
   'param_2': {
      'name': '[]', # Use internal name of the parameter,
      'type': '[]', # 'bool', 'int', 'float' or 'string'
      'val': []  # bool/int/float/string, single-valued or array
            },
   ...      
            
        }
```

- Part 4: Save the results, choose based on the Part2 (Whole_data, Channel_data)
  
**Whole_data**

If there is one output, use the below structure.
```
# Part 4: Whole_data
# Store the result of the function
# You should modify the line codes and replace the placeholders with the internal name of the output
# Change n, the number related to the function that you have defined in main_python file.  
feature['methods'][f'method_{Num_fun}']['outputs'] = {
   'name': [], # Use internal variable name of the output
   'type': [], # output type (string): 'bool', 'int', 'float' or 'string'
   'val': result_n[output]
         }
# Store the warning and success status
# Change n, the number related to the function that you have defined in main_python file.  
feature['methods'][f'method_{Num_fun}']['warning'] = result_n['warning']    
feature['methods'][f'method_{Num_fun}']['success'] = result_n['success'] 
```
If there are multiple outputs, use the below structure.
```
# Store the outputs of the function
# You should modify the line codes and replace the placeholders with the internal name of the output
# Change n, the number related to the function that you have defined in main_python file. 
feature['methods'][f'method_{Num_fun}']['outputs'] = {
   'output_1': {
      'name': [], # Use internal variable name of the output
      'type': [], # output type (string): 'bool', 'int', 'float' or 'string'
      'val': result_n[output_1]
            },
   'output_2': {
      'name': [], # Use internal variable name of the output
      'type': [], # output type (string): 'bool', 'int', 'float' or 'string'
      'val': result_n[output_2]
            },
   ...               
        }
# Store the warning and success status
# Change n, the number related to the function that you have defined in main_python file.  
feature['methods'][f'method_{Num_fun}']['warning'] = result_n['warning']    
feature['methods'][f'method_{Num_fun}']['success'] = result_n['success'] 
```
**Channel_data**

If there is one output, use the below structure.
```
# Part 4: Channel_data
# Store the result of the function
# You should modify the line codes and replace the placeholders with the internal name of the output
feature['methods'][f'method_{Num_fun}']['outputs'] = {
   'name': [], # Use internal variable name of the output
   'type': [], # output type (string): 'bool', 'int', 'float' or 'string'
   'val': output
         }
# Store the warning and success status (Don't change) 
feature['methods'][f'method_{Num_fun}']['warning'] = warning    
feature['methods'][f'method_{Num_fun}']['success'] = success 
```
If there are multiple outputs, use the below structure.
```
# Store the outputs of the function
# You should modify the line codes and replace the placeholders with the internal name of the output
feature['methods'][f'method_{Num_fun}']['outputs'] = {
   'output_1': {
      'name': [], # Use internal variable name of the output
      'type': [], # output type (string): 'bool', 'int', 'float' or 'string'
      'val': output1
            },
   'output_2': {
      'name': [], # Use internal variable name of the output
      'type': [], # output type (string): 'bool', 'int', 'float' or 'string'
      'val': output2
            },
   ...               
        }
# Store the warning and success status
# Change n, the number related to the function that you have defined in main_python file.  
feature['methods'][f'method_{Num_fun}']['warning'] = warning    
feature['methods'][f'method_{Num_fun}']['success'] = success
```
4. For each feature extraction function, copy the four parts and modify them according to the function's name and its outputs.
5. Open the [function_callers](function_callers.py). 
6. Each function has a `functioncaller` file that includes a function to call a specific function. Choose the script based on whether you want to use the Whole_data or Channel_data, and create the corresponding functioncaller_n file where n represents the number related to the function that you have defined in main_python file.
7. Modify `Part 2: Load External Functions` to import all external functions.

**Given the importance of preventing code crashes, it's crucial to have a good familiarity with the input parameters, outputs, and types of them. Therefore, ensure to include both parts of the try and except blocks in the code.**

```Python
# Whole_data
# Change n, the number related to the function that you have defined in main_python file . 
def functioncaller_n(method_name, data, file_name, Params):
    """
    functioncaller_n - Calls a specified function to ... .

    Inputs:
        method_name: Name of the method/function to call
        data: Data on which the function is performed
        file_name: Name of the file (used for error reporting)
        Params: Dictionary containing various parameters, including sampling
                frequency, number of channels, and all input parameters for the feature extraction functions

    Outputs:
        feature: Dictionary indicating the success of the function call and all extracted features
    """

    # Add the line of code "'output': []" based on the number of outputs and replace the placeholders with the internal names of the outputs
    feature = {'success': 0, 'warning': [], 'output1': [], 'output2': [], ...}
    
    try:
        # Call the function 
        function = eval(method_name)
        # Modify Params.var and output based on the function you want to call. Use the internal names of the parameters and outputs. 
        [output1, output2, ...] = function(data, Params.var1, Params.var2, .. )

        # Store results
        # Add the line of code below, based on the number of outputs, and replace the placeholders with the internal names of the outputs
        # Modify the output based on the function you want to call. Use the internal names of the outputs.
        #  ** feature['output'] = output; ** 
        feature['output1'] = output1
        feature['output2'] = output2

        # Indicate success  
        feature['success'] = 1

    except Exception as e:
        # Handle errors during function execution
        feature['warning'] = f"Failed to perform Func: {method_name} for file: {file_name}. Error: {str(e)}"

    return feature
```
```Python
# Channel_data
# k is the index of channel
# Change n, the number related to the function that you have defined in main_pyhton file .  
function channel_feature = functioncaller_n(method_name, data, k, file_name, params)
def functioncaller_n(method_name, data, k, file_name, Params):
    """
    functioncaller_n - Calls a specified function to ...

    Inputs:
        method_name: Name of the method/function to call
        data: Data on which the function is performed
        k: Index of the channel
        file_name: Name of the file (used for error reporting)
        Params: Dictionary containing various parameters, including sampling
                frequency, number of channels, and all input parameters for the feature extraction functions

    Output:
        channel_feature: Dictionary indicating the success of the function call and all extracted features
    """

    # Add the line of code "'output': []" based on the number of outputs and replace the placeholders with the internal names of the outputs
    channel_feature = {'success': 0, 'warning': [], 'output1': [], 'output2': [], ...}

    try:
        # k is the index of channel   
        data_channel = data[k, :]
        # Call the function
        method = eval(method_name)
        # Modify Params.var and output based on the function you want to call. Use the internal names of the parameters and outputs. 
       [output1, output2, ...] = method(data, Params.var1, Params.var2, .. )

        # Indicate success related to each channel
        channel_feature['success'] = 1

        # Store results related to each channel
        # Add the line of code below, based on the number of outputs, and replace the placeholders with the internal names of the outputs
        # Modify the output based on the function you want to call. Use the internal names of the outputs.
        #  ** channel_feature['output'] = output **    
        channel_feature['output1'] = output1
        channel_feature['output2'] = output2
        ... 

    except Exception as e:
        # Handle errors during peak detection
        channel_feature['warning'] = f"Failed to perform Func: {method_name} for file: {file_name}. Error: {str(e)}"
    return channel_feature
```
8. `load` function is essential. Modify `functioncaller_0` to correctly read the data. Each dataset has a specific format, so ensure the loading process is adapted to handle your dataset's format.
9. After modifying the code and creating the relevant `functioncaller_n` files, you can easily run [main_python](main_python.py) in the cluster, AWS, or any other environments.

## Contributors
- [Somayyeh Mousavi](seyedeh.somayyeh.mousavi@emory.edu)
- [Reza Sameni](rsameni@dbmi.emory.edu)