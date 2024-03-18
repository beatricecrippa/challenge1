#include <functional>
#include <vector>
#include <cmath>
#include "helper.hpp"

#ifndef METHOD_HPP
#define METHOD_HPP

class method{
  private:
    
    // struct containing all the elements nneded by the computation
    input _input;
    // contant used to evaluate the finite differences derivative
    static constexpr double h=0.00001;
 
 public:

    //constructor
    method(input &i): _input(i){};
    //constructor by passing as parameter only f and df
    method(functionR  _f, functionRn  _df){
        _input.f=_f;
        _input.df=_df;
    };

    // evaluate gradient in a vector p by finite differences
    vector evaluate_gradient_diff(vector & p){
     vector result;
  
     for(size_t i=0;i<p.size();++i){
        vector v1(p);
        vector v2(p);

        v1[i]+=h;
        v2[i]-=h;
        result.push_back((_input.f(v1)-_input.f(v2))/(2.*h));
     }
    return result;
    }

    // grad computation
     vector compute_grad(vector & xk){
     if (_input.d==Diff::Finite_diff){
        return evaluate_gradient_diff(xk);
     }else if(_input.d==Diff::User_grad){
        return _input.df(xk);
     }
     return evaluate_gradient_diff(xk);
    }
      
    }
    // solving the minimization problem
    vector solve(){

     bool flag=false;
     parameter_type ak=_input.a0;
     unsigned int k=0;
     vector grad;
     vector xold1(_input.start);
     vector xold(_input.start);
     vector xnew(_input.start);

     if(_input.a==Alpha::Armijo && _input.m==Method::Heavy_Ball || _input.a==Alpha::Armijo && _input.m==Method::Nesterov){
      flag=true;
      std::cerr<<"Cannot use Armijo computation for the update of alpha_k and Heavy-Ball or Nesterov method for solving the minimization problem.\n"<<std::endl;
      exit(1);
     }
     while(!flag){
      ++k;
      // xnew computation and gradient evaluation
      grad=compute_grad(xnew);
      xnew=update(xold,xold1,ak,grad,k);

      // check convergence
      flag=check_convergence(xold,xnew,k,grad);
      
      // alpha computation
      ak= compute_alpha(xnew,k);

      //update
      xold=xnew;
      }
     std::cout<<"\nNumber of interations: "<<k<<"\n"<<std::endl;
    return xnew;
  }


    //check convergence
   bool check_convergence(vector & xold,vector & xnew, unsigned int k,vector & grad)const{
        return k>_input.it or euclidean_norm(xold-xnew)<_input.eps_s or euclidean_norm(grad)<_input.eps_r;
    }


    //compute alpha_k - Armijo rule
   parameter_type Armijo(vector & xk){
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