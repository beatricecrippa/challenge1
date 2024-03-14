#include <functional>
#include<vector>
#include<cmath>
#include "helper.hpp"

#ifndef METHOD_HPP
#define METHOD_HPP

class method{
  private:

    input _input;
    static constexpr double h=0.0001;
 
 public:

    //constructor
    method(input i): _input(i){};
    //constructor by passing as parameter only f and sf
    method(functionR _f, functionRn _df){
        _input.f=_f;
        _input.df=_df;
    };
    // evaluate gradient in a vector p by finite differences
    vector evaluate_gradient_diff(vector & p){
     vector result;
     vector grad;
  
     for(size_t i=0;i<p.size();++i){
        vector v1(p);
        vector v2(p);

        v1[i]+=h;
        v2[i]-=h;

        result[i]=(_input.f(v1)-_input.f(v2)/2.*h);
     }
    return result;
    }

    // grad computation
    vector compute_grad(vector & xk){
     if(_input.d==Diff::Finite_diff){
        return evaluate_gradient_diff(xk);
     }else if(_input.d==Diff::User_grad){
        return _input.df(xk);
     }
    return evaluate_gradient_diff(xk);
    }

    // compute the minimizing x
    vector solve(){
     bool flag=false;
     parameter_type ak=_input.a0;
     unsigned int k=0;
     vector grad;
     vector xold(_input.start);
     vector xnew(_input.start);
     while(!flag){

      // alpha computation
      ak= compute_alpha(xnew,k);
    
      // xnew computation and gradient evaluation
      grad=compute_grad(xnew);
      xnew=xold-ak*grad;

      // check convergence
      flag=check_convergence(xold,xnew,k,grad);
      //update
      xold=xnew;
      }
    return xnew;
  }

    

    //check convergence
   bool method::check_convergence(vector & xold,vector & xnew, unsigned int k,vector & grad)const{
        return k>_input.it or euclidean_norm(xold-xnew)<_input.eps_s or euclidean_norm(grad)<_input.eps_r;
    }


    //compute alpha_k - Armijo rule
   parameter_type method::Armijo(vector & xk){
      parameter_type ak=_input.a0;
      value_type gradient_norm_squared=pow(euclidean_norm((evaluate_gradient_diff(xk))),2);
      while(_input.f(xk)-_input.f(xk-ak*evaluate_gradient_diff(xk))<_input.sigma*ak*gradient_norm_squared){
      ak/=2;
    }
    return ak;
    }
    
    //compute alpha_k - exponential decay
    parameter_type exponential_decay(unsigned int k){return _input.a0*exp(-_input.mu*k);}

    //compute alpha_k - inverse decay
    parameter_type inverse_decay(unsigned int k){return _input.a0/(1+_input.mu*k);}

    parameter_type compute_alpha(vector & xk,unsigned int k){
    if (_input.a == Alpha::Armijo) {
        return Armijo(xk);
    } else if (_input.a == Alpha::exponential_decay) {
        return exponential_decay(k);
    } else if (_input.a == Alpha::inverse_decay) {
        return inverse_decay(k);
    }
    return Armijo(xk);

};

#endif
};