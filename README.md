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
| gradient | ![equation](http://www.sciweavers.org/tex2img.php?eq=x_%7Bk%2B1%7D%3Dx_%7Bk%7D-%5Calpha_%7Bk%7D%2A%5Cnabla%20f%28x_%7Bk%7D%29&bc=White&fc=Black&im=jpg&fs=12&ff=mathdesign&edit=0) |
| Heavy-Ball | ![equation](http://www.sciweavers.org/tex2img.php?eq=x_%7Bk%2B1%7D%3Dx_%7Bk%7D-%5Calpha_%7Bk%7D%2A%5Cnabla%20f%28x_%7Bk%7D%29%2B%5Ceta%2A%28x_%7Bk%7D-x_%7Bk-1%7D%29&bc=White&fc=Black&im=jpg&fs=12&ff=mathdesign&edit=0) |
| Nesterov |![equation](http://www.sciweavers.org/tex2img.php?eq=x_%7Bk%2B1%7D%3Dx_%7Bk%7D%2B%5Ceta%2A%28x_%7Bk%7D-x_%7Bk-1%7D%29-%5Calpha_%7Bk%7D%2A%5Cnabla%20f%28x_%7Bk%7D%2B%5Ceta%2A%28x_%7Bk%7D-x_%7Bk-1%7D%29%29&bc=White&fc=Black&im=jpg&fs=12&ff=mathdesign&edit=0) |

Add one of the following keywords  after ```./main``` to decide the method to use to compute the learning rate at each iteration 


| Key word             | meaning                                                                |
| ----------------- | ------------------------------------------------------------------ |
| inverse_decay| ![equation](http://www.sciweavers.org/tex2img.php?eq=%5Calpha_%7Bk%7D%3D%5Cfrac%7B%5Calpha_%7B0%7D%7D%7B%281%20%2B%20%5Cmu%2Ak%29%7D&bc=White&fc=Black&im=jpg&fs=12&ff=mathdesign&edit=0) |
| exponential_decay | ![equation](http://www.sciweavers.org/tex2img.php?eq=%5Calpha_%7Bk%7D%3D%7B%5Calpha_%7B0%7D%7D%2Ae%5E%7B%281%20%2B%20%5Cmu%2Ak%29%7D&bc=White&fc=Black&im=jpg&fs=12&ff=mathdesign&edit=0) |
| Armijo | ![equation](http://www.sciweavers.org/tex2img.php?eq=%5Calpha_%7B0%7D%3D%5Cfrac%7B%5Calpha_%7B0%7D%7D%7B2%7D&bc=White&fc=Black&im=jpg&fs=12&ff=mathdesign&edit=0) until ![equation](http://www.sciweavers.org/tex2img.php?eq=f%28x_%7Bk%7D%29-f%28x_%7Bk%7D-%5Calpha_%7B0%7D%2A%5Cnabla%20f%28x_%7Bk%7D%29%29%20%5Cgeq%20%5Csigma%2A%5Calpha_%7B0%7D%20%5Cparallel%20%5Cnabla%2Af%28x_%7Bk%7D%29%29%20%5Cparallel%20%5E2%0A&bc=White&fc=Black&im=jpg&fs=12&ff=mathdesign&edit=0) is reached|

Add one of the following keywords after ```./main``` to decide the method to compute the gradient  


| Key word             | meaning                                                                |
| ----------------- | ------------------------------------------------------------------ |
| user_grad | the gradient is given by the user |
| finite_diff | the gradient is calculated with the finite differences method |
