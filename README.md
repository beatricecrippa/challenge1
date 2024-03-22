# Challenge 1: A gradient method for the minimization of a multivariate function

### Introduction
Optimization algorithms play a crucial role in numerous computational tasks, ranging from machine learning model training to scientific simulations. This repository hosts an implementation of an optimization algorithm designed to efficiently discover the minimum of a specified function.
### Features
By modifying line 23 in the main file, users can select various configurations, including the method for computing the gradient, the approach for computing the learning rate (alpha_k), and the method used to solve the minimization problem.

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

Using the Getpot library, the user also has the possibility to set the numerical parameters required for solving the minimization problem.

### Repository Contents
+ include: Contains header files, including method.hpp and helper.hpp.
+ src: Contains source code files, with main.cpp being the primary file.
+ Makefile: Provides basic setup for compiling the C++ project.

Contents of helper.hpp

+ Typedefs
+ Enum classes for configuration
+ Struct aggregating computation parameters
+ Function to read from command line and updating the corresponding field in the struct
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
To employ different methods (e.g., Heavy-Ball or Nesterov, use the gradient provided by the user, compute alpha with exponential or inverse decay), modify the solve function template brackets in src/main.cpp at line 23 to specify the desired method, according to the enum classes showed above.

+ First field: Gradient computation method
+ Second field: Learning rate (alpha) computation method
+ Third field: Method for solving the minimization problem


To set the parameters required by the algorithm, it is possible, thanks to the GetPot library, to pass the parameters through the command line, for example:

```bash
    ./main it=10000
```

Otherwise, if only `./main` is executed, the parameters will be fixed by default (for example, by default, the maximum number of iterations is set to it=1000).

Thanks to the following function and the GetPot library, the values passed through the command line will be saved in the corresponding field in the 'input' struct, which contains all the parameters. Below you can see the parameters set by default if not specified in the command line, along with the respective entry to indicate in the command line if you want to modify the value.

```C++

// read struct values by command line
input read_cl(GetPot cl){
    input i;

    // control on the residual
    i.eps_r=cl("eps_r",1e-6);
    // control on the step length
    i.eps_s=cl("eps_s",1e-6);
    // maximum number of iterations
    i.it=cl("it",1000);
    
    // parameters needed by the Armijo rule
    i.sigma=cl("sigma",0.1);
    i.mu=cl("mu",0.2);

    // parameter needed by the Heavy-Ball method
    i.eta=cl("eta",0.7);

    // initial guesses
    i.a0=cl("a0",0.1);

    return i;

}
```



### Lessons Learned
+ Importance of documenting with README.md
+ Utilization of templates and if constexpr statements
+ Utilization of GetPot library                                 
