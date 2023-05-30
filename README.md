# Fully FVM and ALE-based FfSI Solver with a Novel Hybrid Grid
This repository contains the codes for Arbitrary Lagrangian Eulerian (ALE)-based in-house Fluid flexible-Structure Interaction (FfSI) solver with Novel Hybrid AMM-MRR-SEMM Grid Generation Methodology on a Curvilinear Structured Grid with Fully FVM-based discretisation.

# Software requirements
This solver needs:

- gcc
- mpi (for parallel implementation)

# How to install the required packages (on a Linux system)

To install gcc (compiler for C code)

```bash
sudo apt install build-essential
```

To install mpi (OpenMPI)

```bash
sudo apt install mpich
```

# How to complie and run the code

To compile the code

```bash
./compile.sh
```
To run this code

```bash
./output
```
