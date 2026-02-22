# Introduction
This project is created to supplement the following article: (url)
The goal of this project is to demonstrate the effectiveness of reduction of variables in coupling matrix synthesis problems.
Syntheses of extended-box topology coupling matrices of order 8, 14, and 22 are implemented.

# Contents
source code
- general functions: Contains functions required for the synthesis of redundant topologies, and those that implement the similarity transformations.
- extended_box_specific functions: Contains functions specific to the synthesis of extended-box topology coupling matrices, including the implementations of the algorithms present in (url).
- synthesis: Contains the main codes that execute the synthesis procedure.
M_examples: Contains coupling matrices used as examples in (url), together with the free variables obtained from the conducted synthesis procedure.


# How to use
1. Prepare a folded-topology matrix of the desired frequency response, and place it under the folder "source code" as a txt file named "M.txt".
2. Go to ../source code/synthesis/ and open the code "synthesis_through_reduction_C.m" (which implements the synthesis using formulation C in the mentioned article).
3. The settings of the provided code are adjusted for 8th-order case. Change the settings for the desired order. (setting examples for 8th, 14th, and 22nd cases are provided)
4. Run the code, and wait until the solutions are found. The obtained solution matrices are in the variable "sol_M", and the corresponding free variables are in "sol".

<p align="center">
  <img src="M_examples/sword_to_modified_shoelace_30th.gif" alt="Sequence of Givens bijective similarity transformations from 30th-order sword topology to modified-shoelace topology">
</p

*Sequence of Givens bijective similarity transformations from 30th-order sword topology to modified-shoelace topology*



