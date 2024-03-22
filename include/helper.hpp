#include <vector>
#include <cmath>
#include<functional> 
#include <iostream>
#include <map>
#include <iomanip>



// type definitions
typedef double value_type;
typedef std::function<value_type(std::vector<value_type>)> functionR; // R^n -> R function
typedef std::function<std::vector<value_type>(std::vector<value_type>)>  functionRn;  //  R^n -> R^n gradient

#ifndef HELPER_HPP
#define HELPER_HPP

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

// struct aggregating all parameters
struct input{
    
    // function and its gradient, given by the user
    functionR f;
    functionRn df;

    // value of the step size to compute the derivative
    const value_type h=0.00001;

    // control on the residual
    value_type eps_r=1e-6;
    // control on the step length
    value_type eps_s=1e-6;
    // maximum number of iterations
    unsigned int it=1000;
    
    // parameters needed by the Armijo rule
    value_type sigma=0.1;
    value_type mu=0.2;

    // parameter needed by the Heavy-Ball method
    value_type eta=0.7;

    // initial guesses
    value_type a0=0.1;
    std::vector<value_type> start{0.,0.};
};

// print a vector
void print_vector(const std::vector<value_type> & x){
    std::cout<<"[ ";
    for(std::size_t i=0;i<x.size();++i){
        std::cout<<" "<<x[i]<<" ";
    }
    std::cout<<" ]"<<std::endl;

}

// print all elements in the struct containing parameters
template<Diff D,Alpha A,Mode M>
void print_struct(input i){
    std::cout<<"\n\t---Parameters---\n";
    std::cout<<"\nControl on the residual:\t\t\teps_r = "<<i.eps_r;
    std::cout<<"\nControl on the step length:\t\t\teps_s = "<<i.eps_s;
    std::cout<<"\nMaximum number of iterations:\t\t\tk_max = "<<i.it;
    std::cout<<"\nParameters needed by the Armijo rule:\t\tsigma = "<<i.sigma<<"\t\tmu = "<<i.mu;
    std::cout<<"\nParameter needed by the Heavy-Ball method:\teta = "<<i.eta;
    std::cout<<"\nInitial guesses:\t\t\t\talpha_zero = "<<i.a0<<"\tstart_vector = ";
    print_vector(i.start);
    std::cout<<"\n\t---Problem resolution choices---\n";
    std::cout<<"\nMethod for solving the minimization problem: ";
    if constexpr(M==Mode::gradient){
        std::cout<<"\tGradient";
    }else if constexpr(M==Mode::Heavy_Ball){
        std::cout<<"\tHeavy-Ball";
    }else if constexpr(M==Mode::Nesterov){
        std::cout<<"\tNesterov";
    }
    std::cout<<"\nMethod for update alpha: ";
    if constexpr(A==Alpha::Armijo){
        std::cout<<"\t\t\tArmijo";
    }else if constexpr(A==Alpha::exponential_decay){
        std::cout<<"\t\t\texponential decay";
    }else if constexpr(A==Alpha::inverse_decay){
        std::cout<<"\t\t\tinverse decay";
    }
    std::cout<<"\nMethod for computing the gradient: ";
    if constexpr(D==Diff::finite_diff){
        std::cout<<"\t\tFinite differences";
    }else if constexpr(D==Diff::user_grad){
        std::cout<<"\t\tGradient given by the user";
    }
    std::cout<<"\n-------------------------------------------------\n";

}

// compute the Euclidean norm
value_type euclidean_norm(const std::vector<value_type> & v){
    value_type norm=0;
    for(size_t i=0;i<v.size();++i){
        norm+=std::pow(v[i],2.);
    }
    return sqrt(norm);
}

// operator *
std::vector<value_type> operator*(const value_type scalar,const std::vector<value_type> & a){
    std::vector<value_type> result(a.size(),0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i]=(scalar*a[i]);
    }
    return result;
}

// operator -
std::vector<value_type> operator-(const std::vector<value_type> & lhs,const std::vector<value_type> & rhs) {
    std::vector<value_type> result(lhs.size(),0);
    if (lhs.size() != rhs.size()) {
        std::cerr << "\nerror: cannot compute the subtraction between two vector with different size!\n"
                    << std::endl;
        exit(1);
    }
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result[i]=(lhs[i] - rhs[i]);
    }
    return result;
}

// operator +
std::vector<value_type> operator+(const std::vector<value_type> & lhs,const std::vector<value_type> & rhs) {
    std::vector<value_type> result(lhs.size());
    if (lhs.size() != rhs.size()) {
        std::cerr << "\nerror: cannot compute the sum between two vector with different size!\n"
                    << std::endl;
        exit(1);
    }
    for (std::size_t i = 0; i < lhs.size(); ++i) {
        result[i]=(lhs[i] + rhs[i]);
    }
    return result;
}

#endif