"""
Overview
----------
This is the user documenation for the optional Python tools, which are included in the SpinParser software. 
The developer documentation for the SpinParser software itself can be found [here](https://fbuessen.github.io/SpinParser/doc). 
See also the usage examples in the SpinParser [quick start guide](https://github.com/fbuessen/SpinParser#evaluate-spinparser-output-and-measurements). 

Installation
----------
In order to make the Python tools available, simply add the directory `opt/python` in your SpinParser installation to the `$PYTHONPATH` environment variable, such that the scripts can be found by the Python installation.

The following Python packages are prerequisite:

- `h5py`
- `numpy`
- `matplotlib`

The prerequisites can be installed via invoking `python -m pip install h5py numpy matplotlib`.

Once installed, the Python tools provide access to two modules: `spinparser.ldf` and `spinparser.obs`. 
Their functionality is documented below. 
"""