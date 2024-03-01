# Iterative Solutions of Linear Equation Systems (LES)

## Overview
The **Iterative Solutions of LES** project provides MATLAB functions for solving linear equation systems using iterative methods. The implemented methods include Jacobi iteration, Gauss-Seidel iteration, and an examination of the parameter for Smoothed Jacobi iteration. These functions aim to approximate the solution vector of a linear equation system.

- Tests are available within the source codes 

## Functions

### 1. Jacobi Iteration (`jacobi.m`)
Compute the approximation of the solution vector using Jacobi iteration.

**Input:**
- Matrix of the linear equation system A
- Vector for the right side of the equation system b

**Output:**
- Approximation of the solution vector x

### 2. Gauss-Seidel Iteration (`gaussseid.m`)
Compute the approximation of the solution vector using Gauss-Seidel iteration.

**Input:**
- Matrix of the linear equation system A
- Vector for the right side of the equation system b

**Output:**
- Approximation of the solution vector x

### 3. Smoothed Jacobi Iteration (`jomega.m`)
Examine the parameter of Smoothed Jacobi iteration, including drawing the graph of eigenvalues based on the parameter ω.

**Input:**
- Matrix of the linear equation system A

**Output:**
- Spectral radius
- Optimal parameter for convergence
- Interval of convergence
