# Parallelising Conway's Game of Life

This repository allows the user to execute Conway's Game of Life using distributed parallelism with the Message Passing Interface (MPI). This project demonstrates how non-blocking communications between cores can be used to calculate Conway's Game of Life with drastically improved performance over serialized versions of the software. The synopsis of the findings are found in the `xxxx` report.

## Getting Started

This program is predominantly written in C++ with pre-/post-processors written in Python.

It is vital that the user is equipped with a workspace containing MPI, preferable in the Microsoft Visual Studio Community (MVSC) integrated development environment, although help is given below for Linux and iOS.

To get started, simply clone this GitHub repo onto your local machine.

### Prerequisites

C++ and Python.

For Windows (MVSC):
* MVSC Environment pre-equipped with MPI includes and libraries:
  - https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi
  - Download `msmpisdk.mpi` and `msmpisetup.exe`
  - Add location of `mpiexec.exe` to your system PATH
  - Setup  your MVSC environment by:
  -- Adding "C:\Program Files (x86)\Microsoft SDKs\Include" to your workspace Include Directories
  -- Adding “C:\ Program Files (x86)\Microsoft SDKs \Lib\x64" to your workspace Library Directories
  -- Adding "msmpi.lib" to your Linked Input Additional Dependencies.
 
 * iOS/Linux:
    - mpicxx installed on BASH terminal.

## Instructions

Python scripts act as pre- and post-processors to set-up the environment and visualize the results, respectively. Both are written in python scripts, as well as in Jupyter form.

Pre-processor:
* The preprocessor is run to clear the environment for a new run of the game. The user can also input an image (stored in the `./Images/` folder) and run the game with this as their starting image, the user must, however change the file name and `imax_global` and `jmax_global` (domain size) as prompted. In the case of using an image, it is vital that the user must set the `P` parameters (number of cores) and run the script using this number.

Life.cpp:
* The user is required to change the hyper-parameters at the top of the file to play the game they way they wish.
For Windows users: build solution and then call `mpiexec -n P Life.exe` in the `./x64/Release/` directory, where P represents the number of cores to run on.

For iOS/Linux users: compile solution using `mpicxx Life.cpp boundaries.cpp` from within the `./Life/` directory and then call `mpiexec -n P ./a.out`, where P represents the number of cores to run on.

Post-processor:
* The post-processor generates a gif to visualize the results with a movie, along with imaging how the domain has been partitioned. The script runs almost automatically, with the exception that the path must be specified according to your OS (templates provided in the script and the user can delete accordingly). The gif generated is found in the home directory: `game.mp4`.

## Running the tests

To ensure the program runs correctly, the user is encouraged to run the tests to ensure that the games runs correctly, the boundary conditions are correct and the cores are communicating.

Simply set the `test` flag to true at the top of the `Life.cpp` file, where the global environment variables are set, alongside the desired global rows and columns (`imax_global` & `jmax_global`). The user can then compile the solution and run.

For Windows users: build solution and then call `mpiexec -n P Life.exe` in the `./x64/Release/` directory, where P represents the number of cores to run on. For iOS/Linux users: compile solution using `mpicxx Life.cpp boundaries.cpp` from within the `./Life/` directory and then call `mpiexec -n P ./a.out`, where P represents the number of cores to run on.

The post-processor may then be run. The post-processor runs without any required user input except for initially defining the path to the `meta` and `data` directories. (e.g. Windows: `./x64/Release/meta` | `./x64/Release/data` / iOS and Linux: `./Life/meta` | `./Life/data`).

## Built With

* Microsoft Visual Studio Community (MVSC)
* Message Passing Interface (MPI)
* C++
* Python

## Authors

* **Rory Johnston** - [rej19](https://github.com/acse-rej19)

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Stephen Neethling and his team of demonstrators for useful advice and help on how to approach this problem.
