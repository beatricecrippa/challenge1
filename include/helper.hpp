#include <vector>
#include <cmath>
#include <functional>
#include<iostream>


typedef double value_type;
typedef std::vector<value_type> vector;
typedef float parameter_type;

typedef std::function<value_type(vector)> functionR; // R^n -> R function
typedef std::function<vector(vector)>  functionRn;  //  R^n -> R^n gradient

#ifndef HELPER_HPP
#define HELPER_HPP

// this class allows to choose a method for the computation of alpha k
enum class Alpha{
    Armijo,
    exponential_decay,
    inverse_decay,
};


// this class allows to choose the computational method for the gradient
enum class Diff{
    Finite_diff,
    User_grad
};

// this class allow to choose the method for solving the problem of minimization
enum class Method{
    Gradient,
    Heavy_Ball,
    Nesterov
};

struct input{
    
    functionR f;
    functionRn df;

    // control on the residual
    parameter_type eps_r=1e-06;
    // control on the step length
    parameter_type eps_s=1e-06;

    //maximum number of iterations
    unsigned int it=1000;
    
    // parameters needed by the Armijo rule
    parameter_type sigma=0.3;
    parameter_type mu=0.2;

    // parameter needed by the Heavy-Ball method
    parameter_type eta=0.9;

    // initial guess
    parameter_type a0=1;
    vector start{0,0};

    //Alpha computation mode
    Alpha a=Alpha::Armijo;

    //Gradient computation mode
    Diff d=Diff::Finite_diff;
    
    //Method for solving the minimization problem
    Method m=Method::Heavy_Ball;
};

//Euclidean norm
value_type euclidean_norm(vector v){
    value_type norm=0;
    for(size_t i=0;i<v.size();++i){
        norm+=std::pow(v[i],2);
    }
    return sqrt(norm);
}

//operator *
vector operator*(parameter_type scalar,vector a){
    vector result;
    for (std::size_t i = 0; i < a.size(); ++i) {
        result.push_back(scalar*a[i]);
    }
    return result;
}

//operator -
vector operator-(vector lhs,vector rhs) {
    vector result;
    if (lhs.size() != rhs.size()) {
        std::cerr << "\nerror: cannot compute the subtraction between two vector with different size!\n"
                    << std::endl;
        exit(1);
    }
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result.push_back(lhs[i] - rhs[i]);
    }
    return result;
}

    //operator +
vector operator+(vector lhs,vector rhs) {
    vector result;
    if (lhs.size() != rhs.size()) {
        std::cerr << "\nerror: cannot compute the subtraction between two vector with different size!\n"
                    << std::endl;
        exit(1);
    }
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result.push_back(lhs[i] + rhs[i]);
    }
    return result;
}


    //print vector
void print(vector & x){
    std::cout<<"\nThe vector is: [ ";
    for(std::size_t i=0;i<x.size();++i){
        std::cout<<" "<<x[i]<<" ";
    }
    std::cout<<" ]"<<std::endl;

}

#endif
