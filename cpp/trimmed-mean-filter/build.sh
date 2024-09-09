#!/bin/bash

# Compile the source files and generate the executable
g++ -o trimmed_filter trimmed_filter_main.cpp trimmed_filter.cpp list.cpp

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful. Executable 'trimmed_filter' generated."
else
    echo "Compilation failed."
fi