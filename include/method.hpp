#include <functional>
#include <vector>
#include <cmath>
#include "helper.hpp"

#ifndef METHOD_HPP
#define METHOD_HPP

// this class allows to perform the computations to resolve the minimization problem
class method{
  private:
    // struct containing all the parameters needed by the computation
    input _input;
    // constant used to evaluate the finite differences derivative
    static constexpr double h=0.00001;
 
 public:
    // constructor
    method()=default;
    // constructor by passing the struct
    method(input &i): _input(i){print_struct(_input);};

    // solve the problem 
    std::vector<value_type> solve(){
     bool flag=false;
     parameter_type ak=_input.a0;
     unsigned int k=0;
     std::vector<value_type> grad((_input.start).size(),0);
     std::vector<value_type> grad_y((_input.start).size(),0);
     std::vector<value_type> xold(_input.start);
     std::vector<value_type> xnew(_input.start);
     std::vector<value_type> xold1(_input.start);

     if(_input.a==Alpha::Armijo && (_input.m==Mode::Heavy_Ball || _input.m==Mode::Nesterov)){
      flag=true;
      std::cerr<<"Cannot use Armijo computation for the update of alpha_k and Heavy-Ball method for solving the minimization problem.\n"<<std::endl;
      exit(1);
     }
     while(!flag){

    // update number of iterations
    ++k;

    //gradient evaluation
    switch(_input.d){
      case Diff::finite_diff:
        for(size_t i=0;i<xold.size();++i){
        std::vector<value_type> v1(xold),v2(xold);
        v1[i]+=h;
        v2[i]-=h;
        grad[i]=((_input.f(v1)-_input.f(v2))/(2.*h));
        }
      break;
      case Diff::user_grad:
        grad=_input.df(xold);
      break;
       
      }

      // xnew update
      switch(_input.m){
        case Mode::Heavy_Ball:
        if(k==1){
        xnew=xold-ak*grad;
      }else{
        xnew=xold-ak*grad+_input.eta*(xold-xold1);
      }
      break;
      case Mode::gradient:
        xnew=xold-ak*grad;
        break;
      case Mode::Nesterov:
      // gradient y evaluation
      switch(_input.d){
       case Diff::finite_diff:
         for(size_t i=0;i<(xold+_input.eta*(xold-xold1)).size();++i){
          std::vector<value_type> v1(xold+_input.eta*(xold-xold1)),v2(xold+_input.eta*(xold-xold1));
          v1[i]+=h;
          v2[i]-=h;
          grad_y[i]=(_input.f(v1)-_input.f(v2))/(2.*h);
         }
       break;
       case Diff::user_grad:
        grad_y=_input.df(xold+_input.eta*(xold-xold1));
        break;
      }
      
      // alpha computation
      switch(_input.a){
        case Alpha::Armijo:
        ak=_input.a0;
        while(_input.f(xnew)-_input.f(xnew-ak*grad)<_input.sigma*ak*pow(euclidean_norm((grad)),2)){
        ak/=2;
        }
        break;
        case Alpha::exponential_decay:
        ak=_input.a0*exp(-(_input.mu*static_cast<parameter_type>(k)));
        break;
        case Alpha::inverse_decay:
        ak=_input.a0/(1+(_input.mu*static_cast<parameter_type>(k)));
        break;
      }

      // xnew computation
      if(k==1){
        xnew=xold-ak*grad;
      }else{
        xnew=xold+_input.eta*(xold-xold1)-ak*grad_y;
      }
      break;
      }
      

      // check convergence
      flag=check_convergence(xold,xnew,k,grad);
      

      // update
      xold1=xold;
      xold=xnew;

      }
     std::cout<<"Number of interations: "<<k<<"\n"<<std::endl;
    return xnew;
  }
  


  // check convergence
   bool check_convergence(std::vector<value_type> & xold,std::vector<value_type> & xnew, unsigned int k,std::vector<value_type> & grad)const{
        if(k>_input.it){
          std::cout<<"\nLimit of the max number of iterations: no convergence is reached.\n"<<std::endl;
        }else if(euclidean_norm(xold-xnew)<_input.eps_s){
          std::cout<<"\nThe convergence is reached:  ∥ x_k+1 - x_k ∥ < ε_s\n"<<std::endl;
        }else if(euclidean_norm(grad)<_input.eps_r){
          std::cout<<"\nThe convergence is reached: ∥ ∇ f(x_k) ∥ < ε_r\n"<<std::endl;
        }
        return k>_input.it or euclidean_norm(xold-xnew)<_input.eps_s or euclidean_norm(grad)<_input.eps_r;
    }

};

#endif
