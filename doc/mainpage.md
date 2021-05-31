# SpinParser developer documentation

This part of the documentation only contains information which is relevant if you are planning to modify or extend the source code of the SpinParser applications. 

The SpinParser application, without modification, is already suited to solve a broad class of spin models on customizable lattice geometries. 
**The user manual, as well as installation instructions, can be found on the project's [Github page](https://github.com/fbuessen/SpinParser)**. 

Potential reasons to modify the source code could include: 
- implementation of custom, specialized flow equations
- implementation of additional measurements or observables
- implementation of addidional information / log output
- Change of the numerical internal integration routines
- Change of the numerical differential equation solver
- Use of some algorithms / libraries in other projects

Since the SpinParser project aims to make pseudofermion functional renormalization group (pf-FRG) calculations available to a broader community *without* the need to write additional code, no large-scale developer guide for the application and its libraries exists. 
However, on the code level, functions and class interfaces are well documented and a minimal amount of context is provided in class definitions. 

For a general understanding of the code, a good starting point for reading are the SpinParser, FrgCore, and EffectiveAction classes. 
The latter two are abstract definitions, which the SpinParser interfaces when solving pf-FRG flow equations; concrete implementations exist e.g. as SU2FrgCore and SU2EffectiveAction. 
The concrete implementations of the EffectiveAction implement e.g. symmetries in the arguments of vertex functions. 
Measurements are performed by interfacing with the Measurement class. 

The classes Lattice, FrequencyDiscretization, and CutoffDiscretization define the mathematical spaces on which the flow equations are solved. At the beginning of the simulation, based on the user-provided lattice spin model, a Lattice object is generated by the LatticeModelFactory. 

Dynamic load balancing is performed by the LoadManager library, which can be flexibly applied also to new flow equation implementations. 