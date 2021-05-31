![](https://github.com/fbuessen/SpinParser/actions/workflows/build.yml/badge.svg)

# SpinParser

SpinParser is a software platform to perform pseudofermion functional renormalization group (pf-FRG) calculations to solve lattice spin models in quantum magnetism. 

The pf-FRG algorithm has first been proposed in [[Reuther and Wölfle (2010)](http://dx.doi.org/10.1103/PhysRevB.81.144410)] for Heisenberg models on two-dimensional lattices geometries. 

The current implementation of the SpinParser is based on a generalized approach [[Buessen et al. (2019)](http://dx.doi.org/10.1103/PhysRevB.100.125164)] which allows to solve general quantum spin Hamiltonians of the form 

<p align="center"><img src="doc/assets/equation_1.png"></p>

where the sum over lattice sites i and j is defined on arbitrary two- or three-dimensional lattices and the spin operator <img src="doc/assets/equation_2.png" style="vertical-align:-4pt"> resembles the <img src="doc/assets/equation_3.png" style="vertical-align:-3pt"> component of a spin-1/2 moment at lattice site i.

## Overview

Owing to the large class of spin models that can be studied within the pf-FRG approach, the aim of this software implementation is to provide fast, yet flexible numerics for the efficient solution of pf-FRG flow equations. 

To allow for easy and universal use, the spin model and the underlying lattice geometry can be defined as plain-text `.xml` files. 
The SpinParser provides a built-in abstraction layer to identify symmetries in the spin lattice model and exploit them in the subsequent solution of the flow equations. 

The numerical core for the solution of the flow equations themselves is designed to run on massively parallel high-performance computing platforms. 
It utilizes a hybrid OpenMP/MPI parallelization scheme to make efficient use of individual shared-memory computing nodes, while still allowing to scale across multiple computing nodes. 

Due to the algebraic nature of the flow equations, being a large set of coupled integro-differential equations, the computation time for individual contributions to the flow equations can vary. 
To mitigate the impact, SpinParser performs dynamic load balancing across the different computing nodes. 
It is thus in principle also possible to efficiently run the code on heterogeneous computing platforms. 

This software contains three different numerical cores which are optimized for different classes of spin models: 
- A numerical core for spin-1/2 SU(N)-symmetric spin models with an extension for spin-S SU(2) models as put forward in [[Baez and Reuther (2017)](http://dx.doi.org/10.1103/PhysRevB.96.045144), [Buessen et al. (2018)](http://dx.doi.org/10.1103/PhysRevB.97.064415)]
- A numerical core for Kitaev-like spin-1/2 models with diagonal interactions as put forward in [[Reuther et al. (2011)](http://dx.doi.org/10.1103/PhysRevB.84.100406)]
- A numerical core for general spin-1/2 models which may include off-diagonal two-spin interactions as put forward in [[Buessen et al. (2019)](http://dx.doi.org/10.1103/PhysRevB.100.125164)]

## Installation

SpinParser needs to be built from source. The build process uses the [cmake](https://cmake.org) build system. 
SpinParser is relatively easy to compile and only depends on a few libraries. 
The build process is described in a step-by-step guide below, exemplified for Ubuntu 20.04 LTS. 

### Prerequisites

You need to ensure that the following software and libraries are installed on the system: 
* Cmake (version 3.16 or newer)
* Boost (version 1.71.0 or newer)
* HDF5 (version 1.10.4 or newer)
* MPI (optional, recommended)
* Doxygen (optional, required for generating documentation files)

Furthermore, in order for the tests to be evaluated correctly, a python installation is required with the following libraries available:
* numpy
* h5py

To ensure that all these libraries are installed, invoke the following OS specific commands in your terminal. (You might want to drop the `python` part if you already have a python installation.)

*Linux:*

```bash
sudo apt install cmake libboost-all-dev libhdf5-dev libopenmpi-dev doxygen graphviz python
```

```bash
python -m pip install numpy h5py
```

*MacOS (via Homebrew):*

```bash
brew install cmake hdf5 boost doxygen graphviz libomp python
```

```bash
python -m pip install numpy h5py
```

The SpinParser might also successfully build with older library versions. 
Especially when tests are disabled (see below), older versions of boost are also compatible; if the generation of documentation files is disabled (see below), older versions of cmake might work. 

### Download sources

Create and enter a working directory, say 
```bash
mkdir spinParser && cd spinParser
```
Next, clone a copy of the source files
```bash
git clone git@github.com:fbuessen/SpinParser.git spinParserSource
```
which should leave you with a directory `spinParserSource`, which, among other files and subdirectories, contains a file `CMakeLists.txt`. 

### Build from source

Now it is time to compile the sources. It is recommended to build the software in a separate directory, which we will name `build`. 
Furthermore, we want to create a directory `install` to which the final compiled software will be moved. 
Therefore, we create two directories and enter the separate build directory by executing (from within the original working directory `spinParser`)
```bash
mkdir build && mkdir install && cd build
```
We instruct cmake to generate the makefiles for our project by invoking
```bash
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../spinParserSource/
```
The above command includes two command line arguments which tell cmake how to build and install the code. 
The argument `-DCMAKE_BUILD_TYPE=Release` specifies that the generated code is intended to be used for production, i.e. the compiler will be instructed to perform code optimization. 
Alternatively, we could set the value to `Debug` which would instruct the compiler to include debug symbols. 
The second argument, `-DCMAKE_INSTALL_PREFIX=../install` defines where the generated code will ultimately be moved. We intend to move it to the `install` directory that we created in the previous step. 

When executing the command, cmake will attempt to locate all relevant libraries and files. 
On a simple Ubuntu installation, it should usually succeed in doing so. 
If this fails, for example because the libraries are installed in non-standard paths, it is possible to provide cmake with additional hints where to search for the libraries. 
For the boost library, such hint would be similar to the command line argument `-DBOOST_ROOT=/path/to/boost/library`. For the HDF5 library, the hint would be `-DHDF5_ROOT=/path/to/hdf5/library`. 

Furthermore, the SpinParser build environment allows you to specify some additional options: 
* `-DSPINPARSER_BUILD_TESTS=OFF` disables building tests (ON by default).
* `-DSPINPARSER_BUILD_DOCUMENTATION=OFF` disables building the documentation / developer's reference (ON by default).
* `-DSPINPARSER_ENABLE_ASSERTIONS=ON` enables some additional memory boundary and consistency checks. Useful when deriving code or building your own extensions, but slows down the application (OFF by default).
* `-DSPINPARSER_DISABLE_MPI=ON` disables MPI parallelization, which allows code building on systems with no MPI library installed. Can be useful for simplified builds for instrumentation or debugging (OFF by default). 
* `-DSPINPARSER_DISABLE_OMP=ON` disables OpenMP parallelization. Can be useful for simplified builds for instrumentation or debugging (OFF by default). 

Once the `cmake` command has completed, the build files have been generated and we are ready to compile, test, and install the code. 
This is done by calling the sequence of commands (from within the `build` directory)
```bash
make -j6
make test
make install
```
Note that in order to speed up the compilation, we included the argument `-j6`, which instructs to launch 6 processes for parallel compilation. The number should be adjusted to the number of available CPU cores in your system. 

When the compilation is finished, the final software has been installed to the the `install` subdirectory, which we enter by executing
```bash
cd ../install
```
The file structure in this directory should look something like 
```
+ install/
	+ bin/
		- spinParser
	+ doc/
		- index.html
		...
	+ examples/
		- squarelattice-AFM.xml
		...
	+ opt/
		+ mathematica/
			- spinparser.m
		+ python/
			+ spinparser
				- ldf.py
				- obs.py 
	+ res/
		- lattices.xml
		- models.xml
		...
```
The file `install/bin/spinParser` is the main executable. 
Some example definitions of lattices and models are included in `install/res/lattices.xml` and `install/res/models.xml`, respectively. 
The developer documentation / API reference is located in `install/doc/index.html`, which you can open in your browser for reading. 

## Quick start
Performing a calculation with the help of SpinParser consists of three stages:
1. Prepare a task file and define the microscopic lattice and spin model.
2. Ensure that your implementation is correct. 
3. Run the calculation. 
4. Evaluate SpinParser output and measurements. 


### Prepare a task file
The task file is an XML document which holds all input required to perform an FRG calculation.
This includes a reference to the underlying lattice, the spin model, information about frequency and cutoff discretizations, and instructions on the measurements to perform. 
An example task file is included in the SpinParser installation under `install/examples/squarelattice-AFM`: 
```XML
<task>
	<parameters>
		<frequency discretization="exponential">
			<min>0.005</min>
			<max>50.0</max>
			<count>32</count>
		</frequency>
		<cutoff discretization="exponential">
			<min>0.1</min>
			<max>50</max>
			<step>0.95</step>
		</cutoff>
		<lattice name="square" range="4"/>
		<model name="square-heisenberg" symmetry="SU2">
			<j>1.0</j>
		</model>
	</parameters>
	<measurements>
		<measurement name="correlation"/>
	</measurements>
</task>
```
This example defines an exponential frequency discretization, specified in the block `<frequency discretization="exponential">`. The distribution is generated symmetrically around zero with 32 frequencies in the range from 0.005 to 50.0, i.e., a total of 64 frequencies are used in the computation. 
Alternatively, an explicit list of values can specified by choosing `discretization="manual"` and providing the values as child nodes `<value>0.005</value> [...] <value>50.0</value>`. 

The cutoff discretization is automatically generated as an exponential distribution <img src="doc/assets/equation_4.png" style="vertical-align:-3pt"> down to the smallest cutoff value <img src="doc/assets/equation_5.png" style="vertical-align:-3pt">, according to the specification in the node `<cutoff discretization="exponential">`. 
Just like in the specification of the frequency discretization, it is also possible to specify `discretization="manual"`.

The lattice graph `<lattice name="square" range="4"/>` will be generated to include all lattice sites up to a four lattice-bond radius around a reference site. The name of the lattice, `square`, is a reference to a lattice definition found elsewhere. The actual lattice definition is found in the resource file `install/res/lattices.xml` file: 
```XML
<unitcell name="square">
	<primitive x="1" y="0" z="0" />
	<primitive x="0" y="1" z="0" />
	<primitive x="0" y="0" z="1" />

	<site x="0" y="0" z="0" />

	<bond from="0" to="0" dx="1" dy="0" dz="0" />
	<bond from="0" to="0" dx="0" dy="1" dz="0" />
</unitcell>
```
The software will scan all `*.xml` files in the directory `../res/` relative to the SpinParser executable for lattice implementations. 
If the directory does not exist, it will scan the directory `../../res/` or the directory where the executable is located. 
Alternatively, the resource search path can be set manually by specifying the `--resourcePath` command line argument when launching SpinParser. 

A lattice definition consists of three `primitive` lattice vectors spanning the unit cell; each is defined by their x, y, and z component. 
The lattice unit cell, in this example, contains one basis site at the origin. Multiple basis sites are in principle possible, in which case they are enumerated by unique IDs according to their order, starting at 0. 
Finally, the lattice graph is generated according to the lattice bonds: each bond connects two lattice sites, `from` and `to` (referenced by their ID), which may either lie within the same unit cell, or be offset by dx, dy or dz unit cells into the direction of the first, second or third lattice vector, respectively. 

Similarly, the spin model `<model name="square-heisenberg" symmetry="SU2">` references the model `square-heisenberg` defined in the file `install/res/models.xml`:
```XML
<model name="square-heisenberg">
	<interaction parameter="j" from="0,0,0,0" to="1,0,0,0" type="heisenberg" />
	<interaction parameter="j" from="0,0,0,0" to="0,1,0,0" type="heisenberg" />
</model>
```
The software will scan all `*.xml` files in the directory `../res/` relative to the SpinParser executable for spin model implementations. 
If the directory does not exist, it will scan the directory `../../res/` or the directory where the executable is located. 
Alternatively, the resource search path can be set manually by specifying the `--resourcePath` command line argument when launching SpinParser. 

The model definition comprises a list two-spin interactions. All interactions for one lattice unit cell need to be specified, the rest is inferred by periodicity of the lattice. 
The interaction is between two lattice sites `from` and `to`, each referenced by a tuple (a1,a2,a3,b), corresponding to the lattice site in unit cell (a1,a2,a3) (in units of the lattice vectors) and basis site ID b. 
The two-spin interaction type in this example is a Heisenberg interaction. 
Any interaction can be specified explicitly by replacing `heisenberg` with `xy` which would correspond to the two-spin interaction <img src="doc/assets/equation_6.png" style="vertical-align:-4pt">. 
The `parameter` name is referenced in the task file to assign a numerical value to the coupling. 
In our example above, in the line `<j>1.0</j>`, it is set to 1.0, with the sign convention such that the interaction is antiferromagnetic. 

The `symmetry` attribute in the model reference of the task file specifies which numerical backend to use. Possible options are `SU2` (compatible with SU(2)-symmetric Heisenberg interactions), `XYZ` (compatible with diagonal interactions) or `TRI` (compatible also with off-diagonal interactions). 
You should generally use the numerical backend with the largest possible symmetry, as this will greatly reduce computation time. 

Finally, the line `<measurement name="correlation"/>` specifies that two-spin correlation measurements should be recorded. 

### Verify the model implementation
To ensure that all interactions have been specified correctly, you can invoke the SpinParser (see also next section) with the command line argument `--debugLattice`, 
```bash
./install/bin/spinParser --debugLattice install/examples/squarelattice-AFM.xml
```
which will not run the actual calculation, but only produce an output file `install/examples/squarelattice-AFM.ldf`, which is an xml-type file which describes all relevant lattice sites with real space coordinates and all lattice bonds and interactions with the 3x3 interaction matrices as labels. 

### Run the calculation
Once the task file is prepared, you can launch the calculation by invoking
```bash
./install/bin/spinParser install/examples/squarelattice-AFM.xml
```
The computation will use OpenMP to utilize the maximum number of available CPU cores according to the environment variable `OMP_NUM_THREADS`. 

If you are running the calculation on a distributed memory machine, you may launch SpinParser in MPI mode, e.g. 
```bash
mpirun -n 8 install/bin/spinParser install/examples/squarelattice-AFM.xml
```
The above command would launch the calculation in a hybrid OpenMP/MPI mode across 8 nodes, using the maximum number of available OpenMP threads on each node. 

As the calculation progresses, an output file `install/examples/squarelattice-AFM.obs` is generated which contains the measurement results as specified in the task file. 

### Evaluate SpinParser output and measurements
The result file `install/examples/squarelattice-AFM.obs` is an hdf5 file which contains the two-spin correlation measurements. 
It contains datasets like `/SU2CorZZ/data/measurement_0/data`, which is a list of two-spin correlations <img src="doc/assets/equation_7.png" style="vertical-align:-4pt"> with lattice sites n in the same order as listed in `/SU2CorZZ/meta/sites`. 
Every dataset is performed at the cutoff value as specified in the attribute `/SU2CorZZ/data/measurement_0/cutoff`. 

The data is now ready to be extracted and analyzed. 

## Developer documentation

The SpinParser application, without modification, is already suited to solve a broad class of spin models on customizable lattice geometries. 
However, some users might wish to modify or extend the code; for this purpose, a developer documentation of the underlying code exists and can be found [here](https://fbuessen.github.io/SpinParser).