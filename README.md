# Challenge 1: A gradient method for the minimization of a multivariate function

### Introduction
Optimization algorithms play a crucial role in numerous computational tasks, ranging from machine learning model training to scientific simulations. This repository hosts an implementation of an optimization algorithm designed to efficiently discover the minimum of a specified function.
### Features
By modifying line 20 in the main file, users can select various configurations, including the method for computing the gradient, the approach for computing the learning rate (alpha_k), and the method used to solve the minimization problem.

With this flexibility, users can adapt the optimization algorithm to suit their specific computational requirements, optimizing performance while minimizing computational resources.

This customization is facilitated by the solve function, which utilizes three enum classes as templates:

```C++
// Enum class enabling the selection of a method for computing alpha_k
enum class Alpha{
    Armijo,
    exponential_decay,
    inverse_decay,
};

// Enum class allowing the choice of gradient computation method

enum class Diff{
    finite_diff,
    user_grad
};

// Enum class enabling the selection of a method for solving the minimization problem
enum class Mode{
    gradient,
    Heavy_Ball,
    Nesterov
};
```
### Repository Contents
+ include: Contains header files, including method.hpp and helper.hpp.
+ src: Contains source code files, with main.cpp being the primary file.
+ Makefile: Provides basic setup for compiling the C++ project.

Contents of helper.hpp

+ Typedefs
+ Enum classes for configuration
+ Struct aggregating computation parameters
+ Functions for vector and struct printing
+ Euclidean norm computation function
+ Operators for vector computations (scalar product, subtraction, and addition)

Contents of method.hpp
+ Function to check convergence of the method
+ Function to compute the gradient
+ Function to solve the minimization problem
### Requirements
+ C++ compiler compatible with C++20 standard
+ Git for cloning the repository
+ Make utility for compiling the code
### Running Locally
Clone the project:
```bash
    git clone git@github.com:irene-fagnani/challenge1.git
```
Compile the code:
```bash
    make
```
Run the challenge:
```bash
    ./main
```
To employ different methods (e.g., Heavy-Ball or Nesterov, use the gradient provided by the user, compute alpha with exponential or inverse decay), modify the solve function template brackets in src/main.cpp at line 20 to specify the desired method, according to the enum classes showed above.

+ First field: Gradient computation method
+ Second field: Learning rate (alpha) computation method
+ Third field: Method for solving the minimization problem
### Lessons Learned
+ Importance of documenting with README.md
+ Utilization of templates and if constexpr statements                                  
