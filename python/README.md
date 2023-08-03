# OSET Tools in Python

## Installing

To use the OSET Python package, follow these steps according to your system's Python release:

### Setting up a Virtual Environment using Conda

1. Open a terminal or command prompt.

2. Create a new virtual environment using Conda. Replace `<venv_name>` with your desired name for the virtual
   environment:

   ```bash
   conda create -n <venv_name> python=3.10 --file requirements.txt
   ```

3. Activate the virtual environment:

   ```bash
   conda activate <venv_name>
   ```

### Setting up a Virtual Environment using Native Python (venv)

1. Open a terminal or command prompt.

2. Navigate to the root directory of your project (where `setup.py` and `requirements.txt` are located).

3. Create a new virtual environment using the built-in `venv` module. Replace `<venv_name>` with your desired name for
   the virtual environment:

   ```bash
   python -m venv <venv_name>
   ```

4. Activate the virtual environment:

    - On Windows:

      ```bash
      <venv_name>\Scripts\activate
      ```

    - On macOS and Linux:

      ```bash
      source <venv_name>/bin/activate
      ```

### Building and installing the OSET package on the Virtual Environment (tested on MacOS)

After activating the venv, make sure to have `poetry` installed. If not, run:

```bash
pip install poetry
```

Now build the OSET package:

```bash
poetry build
```

With the virtual environment activated (as explained above), you can now install the OSET package using `pip`:

```bash
pip install .
```

After installation, you can check if the package is correctly installed in the virtual environment:

```bash
pip list
```

### Building and installing the OSET package on the Virtual Environment (tested on Windows)

Make sure all the package subfolders have a `__init__.py` file. Now build the OSET package:

```bash
python -m build
```

With the virtual environment activated (as explained above), you can now install the OSET package using `pip`:

```bash
pip install .
```

After installation, you can check if the package is correctly installed in the virtual environment:

```bash
pip list
```

### Uninstalling

To uninstall the OSET package, use the following command:

```bash
pip uninstall oset
```

Replace `<venv_name>` with the name of the virtual environment where you installed OSET.

After uninstallation, verify that the package is fully removed by rechecking the installed packages using `pip list`. If
you previously had a different version installed, uninstalling will revert to the previous version.

The virtual environment can be deactivated as follows:

```bash
conda deactivate
```

And it can be deleted using:

```bash
conda remove --name <venv_name> --all
```

## Jupyter

Install Jupyter on the virtual environment:

```bash
conda install jupyterlab
```

After activating the created virtual environment, install its kernel:

```bash
ipython kernel install --user --name=<venv_name>
```

Launch `jupyter-notebook` or `jupyter-lab` after activating the virtual environment. Select the `<venv_name>` kernel
from the Kernels tab.

Check in the Jupyter terminal to make sure that `oset` is installed (`pip list` or `conda list`). If not, install OSET
from the Jupyter terminal.

## Usage

After successfully installing the OSET package in your virtual environment, you can use it in your Python scripts or
Jupyter Notebooks:

```python
from oset.generic.tanh_saturation import tanh_saturation

# Your code using the oset package
```

Please note that you need to activate the virtual environment (as explained above) whenever you want to work with the
OSET package.

For more details and the latest updates, visit the [OSET GitHub repository](https://github.com/alphanumericslab/OSET).

Reza Sameni, 2023  
The Open-Source Electrophysiological Toolbox

# Contribution

## Formatting

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

We use Black for the python part of the repository. Please make sure to format your code with Black before committing.

To install Black, run:

```bash
pip install black black[jupyter]
```

To format your code, run the following command from within the [python](../python) directory of the repository:

```bash
python -m black .
```

This will reformat all the python and jupyter notebook files in the python directory.