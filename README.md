# A Recursive Domain- and Objective-Adaptive Frank-Wolfe Algorithm
This repository contains the Julia code used to generate the evaluation data for the numerical examples in the paper `A Recursive Domain- and Objective-Adaptive Frank-Wolfe Algorithm (May 2026 - Marcel Kaiser, Tobias Sutter)`.

## Repository Structure
The repository contains two examples:

- `Simple_Example/`: first numerical example
- `Linear_Quadratic_Gaussian_Contol/`: second numerical example

Each example can be run independently by executing the corresponding `main.jl` file.

## Requirements

The code is entirely written in Julia.

Some examples require the following Julia packages:

- CairoMakie
- Statistics
- Random
- LinearAlgebra
- Printf
- LaTeXStrings
- JuMP
- MosekTools
- Distributions
- ForwardDiff

The required packages are installed automatically when running the scripts.

The second example uses Mosek through `MosekTools.jl`. Therefore, a working Mosek installation and license are be required.

## Installation
Schritte zum Aktivieren des Environments und Installieren der Pakete.

## How to run

Clone or download this repository.

To run the first example:

```bash
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
