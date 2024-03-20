# Challenge 1
This repository contains the implementation of an optimization algorithm aiming to find the minimum of a given function. 

The user is free to choose:
+ The method used to compute the gradient (by finite differences or given the gradient of the function)
+ The method used to compute the learning rate alpha_k (by inverse decay, exponential decay or the Armijo rule)
+ The method used to solve the minimization problem (by the gradient method, by the Heavy-Ball method or by the Nesterov method). If Nesterov or Heavy-Ball are choosen, the program won't resolve the problem if the computational method for the learning rate is the Armijo rule, because the update direction is not guaranteeded to be a descent direction.

### Requirements
+ C++ compiler compatible with C++20 standard.
+ Git for cloning the repository.
+ Make utility for compiling the code.


## Run Locally

Clone the project

```bash
  git clone git@github.com:irene-fagnani/challenge1.git
```

Compile the code

```bash
  make
```

Run the challenge in the default mode (so applying the Gradient method, the Armijo rule to compute alpha and the finite differences method to compute the gradient)

```bash
  ./main
```
Run other methods (eg. Heavy-Ball or Nesterov, use the gradient given by the user, compute alpha with the exponential or inverse decay rule)

```bash
  ./main Nesterov user_grad inverse_decay
```
Add one of the following keywords after ```./main``` to decide the method to use to solve the minimization problem  


| Key word             | meaning                                                                |
| ----------------- | ------------------------------------------------------------------ |
| gradient | Apply the gradient method |
| Heavy-Ball | Apply the Heavy-Ball method | 
| Nesterov | Apply the Nesterov method|

Add one of the following keywords  after ```./main``` to decide the method to use to compute the learning rate at each iteration 


| Key word             | meaning                                                                |
| ----------------- | ------------------------------------------------------------------ |
| inverse_decay| computing learning rate with inverse decay formula |
| exponential_decay | computing learning rate with exponential decay formula |
| Armijo | computing learning rate with aRmijo rule|

Add one of the following keywords after ```./main``` to decide the method to compute the gradient  


| Key word             | meaning                                                                |
| ----------------- | ------------------------------------------------------------------ |
| user_grad | the gradient is given by the user |
| finite_diff | the gradient is calculated with the finite differences method |

