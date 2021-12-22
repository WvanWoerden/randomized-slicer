# Readme

This repository contains the supplementary material of the following paper: 
- L. Ducas, T. Laarhoven, W. van Woerden, The randomized slicer for CVPP: sharper, faster, smaller, batchier. PKC 2020.

# Contents

- **short_path.cpp:** The C++ code to numerically approximate the optimal reduction path.
```
g++ short_path.cpp -O3 -o short_path -fopenmp
mkdir data
./short_path
```
- **symbolic_verif.ipynb:** A sage notebook verifying the proofs using symbolic algebras.

# Acknowledgments

Leo Ducas was supported by the European Union H2020
Research and Innovation Program Grant 780701 (PROMETHEUS) and the Veni
Innovational Research Grant from NWO under project number 639.021.645.
Thijs Laarhoven was supported by a Veni Innovational Research Grant from
NWO under project number 016.Veni.192.005. Wessel van Woerden was supported by the ERC Advanced Grant 740972 (ALGSTRONGCRYPTO).
