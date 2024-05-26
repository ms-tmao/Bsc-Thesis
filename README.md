# Finding Closest Normal Matrices

This repository contains a MATLAB implementation of algorithms based on A. Ruhe’s paper ‘Closest Normal Matrix Finally Found’ (1987). It is part of a Bachelor's thesis at TU Berlin.

## Routines

The repository contains the following functions:

- closest_normal.m : the main implementation of Ruhe's algorithm
- closest_normal2.m : a 2nd implementation of Ruhe's algorithm using a different pivot strategy
- closest_normal_wsweep.m : w sweeps of the main implementation
- closest_normal2_witer.m : w iterations of the 2nd implementation
- determine_pivot.m : calculates pivots for closest_normal2.m and closest_normal2_witer.m
- obesity_Y.m : constructs an obesity matrix as defined in chapter 2 of the thesis
- spread.m : calculates the spread as in chapter 2 of the thesis
- create_Jordan_matrix.m : constructs Jordan matrices
- nu.m : calculates measures of nonnormality

For the purpose of demonstration, it furthermore contains some example calculations:

- example_Jordan_matrix.m : example featuring a Jordan matrix
- example_Ruhe.m : examples taken from Ruhe's original paper
- example_thesis.m : further examples that can be found in the thesis

## Calculating Closest Normals

To simply calculate a candidate for a closest normal of A, define accuracy d and write 

```matlab
N = closest_normal(A, d)
```

or

```matlab
[N, D, A_new, R, sweep] = closest_normal(A, d)
```

The above array contains
- N the putative closest normal of A
- A_new the Delta H-Matrix that A is transformed into
- D the diagonal of A_new 
- R the unitary matrix transforming A into A_new
- sweep the number of sweeps needed by the algorithm

## Verifying Closest Normals

To use the spread criterion, first calculate A_new (by closest_normal.m, see above) and obesity matrix Y of A_new, and then calculate its spread

```matlab
[N, D, A_new, R, sweep] = closest_normal(A, d);
spread(obesity_Y(A_new))
```
N is a closest normal if the above returns a value smaller than 1.

## Credits

Tong Mao
Last updated: 22 May 2024