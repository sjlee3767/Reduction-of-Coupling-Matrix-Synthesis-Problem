# Introduction
This project is created to supplement the following article: (url)
The goal of this project is to demonstrate the effectiveness of the formulations proposed in the article for solving the coupling matrix synthesis problem.
Syntheses of extended-box topology coupling matrices of order 8, 14, and 22, and the pedagogical 8th-order topology in the appendix are implemented.

# Contents
source code
- general functions: Contains functions required for the synthesis of redundant topologies, and those that implement the similarity transformations.
- topology-specific functions: Contains functions specific to the synthesis of example topologies in the article, including the implementations of the algorithms present in (url).
- synthesis: Contains the main codes that execute the synthesis procedure.
M_examples: Contains coupling matrices used as examples in (url).


# How to use
0. Open Matlab, and add the project folder and all the subfolders to the Matlab path.
1. Prepare a folded-topology matrix of the desired frequency response, and place it under the folder "source code" as a txt file named "M.txt". (The folded-topology matrices used in the article can be found in M_examples/folded.)
2. Go to ../source code/synthesis/ and open the code "eb_synthesis_formulation_E.m" (which implements the synthesis of extended-box topology using formulation E in the mentioned article).
3. Run the code, and wait until the solutions are found. The obtained solution matrices are in the variable "sol_M", and the corresponding free variables are in "sol".

<p align="center">
  <img src="M_examples/sword_to_modified_shoelace_30th.gif"
       alt="Sequence of Givens bijective similarity transformations from 30th-order sword topology to modified-shoelace topology">
  <br>
  <em>Sequence of Givens bijective similarity transformations from 30th-order sword topology to generalized-box topology</em>
</p>

# Notes
- The current synthesis code supports 8th, 14th, and 22nd-order extended-box topology. The synthesis of arbitrary-order generalized-box topology is implemented in generalized_box_nth.

