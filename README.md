# Bayesian Analysis Code for SLI 

## Introduction
This code base is complementary to our End2End shaken lattice project. The C++ component uses shaking sequences that were previously learned using RL to simulate the shaken lattice and record the final momentum space probability distributions $P(p|a,V_0)$, where $p$ is momentum and $a,V_0$ are the acceleration and lattice depths for the system specifically. The Python notebooks are then used to do Bayesian and statistical analysis.

General Notes:

- The shaking sequences are currently hard-coded as constant `std::vectors` into the source files. While not programmatically elegant, this is partly because for the purposes of this project, the shaking sequences should never be changed, even in error, and there are only three sequences.
- The shaking sequences are generated via RL code, which for privacy purposes isn't public yet. It however uses a version of the `Environment` and `Wavepacket` classes shown here.
- The Bayesian updating, Jenson-Shannon divergence calculation and other statistical analyses are all in the Python notebooks. The Python code reads in the files generated by the C++ simulations.
- Despite Github marking this as solely a Jupyter project, the core simulations and physics of the shaken lattice are encapsulated in the C++ code. 

## Structure of the code
The code is structured in the following folders:

- **src** contains the source main files
    - _Compare_ was used to compare the accuracy of our probability distributions to those from the Python code.
    - _Generate_ is used to generate the final momentum probability files for different $(a,V_0)$ values. It does so by directly evolving the `Wavepacket` class.
    - _Simulate\_Momentum_ uses the shaking sequences to evolve the `Momentum` class, thus recording the evolution of different observables with time. (This file is in progress).
- **include** contains the class files: `Environment` and `Wavepacket`.
- **Python** includes the Jupyter notebooks used for generating plots for Bayesian updating as well as assorted plots for the paper.
- **Output** is where the files generated by _Generate_ and _Simulate\_Momentum_ are stored as the code runs.
- **Data** contains past runs from *Output*.
- **build** is for using `CMake` to build the C++ program, and contains binaries. The `CMakeLists.txt` is modified whenever we pick a different target and **src/** file.

