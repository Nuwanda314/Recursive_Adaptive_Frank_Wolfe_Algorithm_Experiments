# A Recursive Domain- and Objective-Adaptive Frank-Wolfe Algorithm
This repository contains the Julia code used to generate the evaluation data for the numerical examples in the paper *A Recursive Domain- and Objective-Adaptive Frank-Wolfe Algorithm (May 2026 - Marcel Kaiser, Tobias Sutter)*.

## Repository Structure
The repository contains two examples:

- `Simple_Example/`
- `Linear_Quadratic_Gaussian_Contol/`

Each example can be run independently by executing the corresponding `main.jl` file in the `src/` folder.

## Requirements

The code is entirely written in Julia. 

The examples require the following packages:

- Simple Example
  - `CairoMakie`
  - `Random`
  - `Statistics`

- Robust Linear Quadratic Gaussian Control
  - `CairoMakie`
  - `Distributions`
  - `ForwardDiff`
  - `JuMP`
  - `LinearAlgebra`
  - `MosekTools`
  - `Printf`
  - `Random`
  - `Statistics`

The required packages are installed automatically when running the corresponding `main.jl` scripts.

The second example uses Mosek through `MosekTools.jl`. Therefore, a working Mosek installation and license are required.

## How to run

Clone or download this repository.

To run the first example:

```
bash
cd example_1
julia Main.jl
```

## Reproducibility
Erklären:
- welche Datei welchen Teil des Papers reproduziert
- wo Parameter eingestellt werden
- wo Outputs landen

## Results
Optional 1–2 Abbildungen mit kurzer Erklärung.

## Citation
BibTeX oder Paper-Zitat.

## Notes
Optional: work in progress, planned extensions, numerical caveats.
