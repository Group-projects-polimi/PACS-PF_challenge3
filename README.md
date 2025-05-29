# PACS challenge 3: parallal Laplace equation solver

- [Functionality](#functionality)
- [Getting started](#getting-started)
- [Implementation](#implementation)
  - [Matrix-free memory management](#matrix-free-memory-management)
  - [Parallel solver restrictions](#parallel-solver-restrictions) 
- [Performance](#performance)

## Functionality

The code implements `serialSolver` and `parallelSolver` classes that numerically approximate the solution to a Laplace equation on a square domain with arbitrary choices of forcing term, Dirichlet boundary conditions and number of grid points, as well as facilities for evaluating the error, exporting the results in _.vtk_ format for visualization in [Paraview][] and automatically generating profiling plots with [Gnuplot][]. All of it comes with a **very fast implementation** that harnesses the joint potential of [OpenMP][] and [MPI][].

- _make libs_ automatically downloads the [Nlohmann JSON][] library and the [Muparser][] library from GitHub, then builds the second one into a shared library and properly links it so that other compilation commands may find it.
- _make_ builds the program with its default settings. Running it with _mpirun -np 4 executable_ will solve the Laplace equation using the selected number of processors according to the data specified inside _data.json_ and export the results in _.vtk_ format.
- _make test_ builds the pogram and performs tests on a model problem with varying number of grid nodes, MPI processes and OpenMp threads, then uses Gnuplot to output the resulting graphs in the _plots_ folder.
- _make debug_ builds the program while enabling many assertions throughout the code that, while disabled by default for efficiency concerns, make it safer to run; indeed if something is not working properly try compiling with this option to see if there's an error in the input or in the sequence of operations or if the code is actually broken.
- _make clean_ and _make doc_ do what they claim.

> ⚠️ **Warning**: Outside of an installation of OpenMP and MPI and optional tools like Paraview and Gnuplot, this project needs the **Nlohmann JSON** library and a shared **Muparser** library to work. _make libs_ downloads their source code, compiles the second into a shared library and puts everything in its proper place. _make libs_ is written to be portable but if it's not working the user is invited to download the dependencies, compile Muparser into a shared library and adhere to the directory structure below.

```bash
laplace_solver/
├── Makefile                                                  # Compilation commands
├── data.json                                                 # Problem data
├── libmuparser.so -> ./muparser/build/libmuparser.so.2.3.5   # Symbolic link
├── libmuparser.so.2 -> ./muparser/build/libmuparser.so.2.3.5 # Symbolic link
├── .include/                                                 # Header files
├── src/                                                      # Source files
├── test/                                                     # Tests' header and source files
├── json/                                                     # Nlohmann JSON parsing library
│   └── single_include/
│       └── nlohmann
└── muparser/                                                 # Muparser library
    └── build/
        └── libmuparser.so.2.3.5
```

## Getting started
To solve a laplace equation just write your _mydata.json_ file or modify the provided one, here it's possible to specify number of grid points, boundary conditions and forcing term among other options. Such files can be passed to the program with _mpirun -np 4 executable mydata.json_ or _./executable mydata.json_.

If instead one wants to modify the source code to obtain a more complex behavior than the default one, then just know that the most convenient way to initialize a `serialSolver` or a `parallelSolver` object is to simply pass a `std::string` containing the name of the file to read data from, such file should adhere to the strcture of _data.json_. Then one can use the solver to approximate the solution to a Laplace problem, print the genrated mesh, print the error (assuming an exact solution was provided) or export the results in _.vtk_ format. A simple example is provided below.

```cpp
parallelSolver s(filename);
s.print_mesh();
s.solve();
s.export_vtk();
```

## Implementation
The code is **doxygen-documented**, the doxygen documentation is the best place to learn more about the inner workings of the code before diving in the source files. What is useful to bring to the reader's attention from the beginning are a couple of details including how the solution matrix is managed in memory and some restrictions of the current implementation. 

### Matrix-free memory management
When executing with more than 1 MPI process the memory load is evenly distributed. Never in the code does a `parallelSolver` hold the whole solution matrix: not even while printing or exporting results. This **matrix-free** and **disributed** approach allows handling very large matrices.

### Parallel solver restrictions
It's not possible to run one single parallel process and ask it to solve a problem with a `parallelSolver` object: the program will crash if compiled with _make_ or an assertion will fail if compiled with _debug_. This is because `parallelSolver`s are made to communicate with each others, if one wishes to have just one process than he should define a `serialSolver`.
It's not possible to run a certain number of parallel processes `N` and solve a problem where the number of grid points is strictly less than `N`: the program will crash if compiled with _make_ or an assertion will fail if compiled with _debug_.

## Performance
To obtain decent performance results it's **vitally important** to set a reasonable number of OpenMp threads and to couple them with a reasonable number of MPI processes, where of course the definition of reasonable depends on one' specific system. To set the nuber of OpenMP threads on Unix machines, simply run `export OMP_NUM_THREADS=4`.
Some results are shown for a low-end machine with these characteristics: intel core i7-8665U CPU @ 1.90GHz, treads per core: 2, cores per socket: 1, socket: 1. Note that these same plots can be obtained by running _make test_.

![time_serial](https://github.com/user-attachments/assets/ed8f9d60-fc7f-4282-94a7-ff6097896130)
![time_parallel](https://github.com/user-attachments/assets/852ff97e-c290-4262-9fc2-5d51a9159e2b)

[Paraview]: https://www.paraview.org/
[Gnuplot]: https://en.wikipedia.org/wiki/Gnuplot
[OpenMP]: https://www.openmp.org/
[MPI]: https://en.wikipedia.org/wiki/Message_Passing_Interface
[Muparser]: https://beltoforion.de/en/muparser/
[Nlohmann JSON]: https://github.com/nlohmann/json
