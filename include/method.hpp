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
    
    void set_eta(parameter_type value){
      _input.eta=value;
    }
    // evaluate gradient in a vector p by finite differences
    std::vector<value_type> evaluate_gradient_diff(const std::vector<value_type> & p){
     std::vector<value_type> result;
     for(size_t i=0;i<p.size();++i){
        std::vector<value_type> v1(p),v2(p);
        v1[i]+=h;
        v2[i]-=h;
        result.push_back((_input.f(v1)-_input.f(v2))/(2.*h));
     }
    return result;
    }

    // grad computation
     std::vector<value_type> compute_grad(const std::vector<value_type> & xk){
     if (_input.d==Diff::Finite_diff){
        return evaluate_gradient_diff(xk);
     }else if(_input.d==Diff::User_grad){
        return _input.df(xk);
     }
     return evaluate_gradient_diff(xk);
    }
    
    // solving the minimization problem
    std::vector<value_type> solve(){
      if(_input.m==Mode::Gradient){
        return gradient_descendant();
      }else if(_input.m==Mode::Heavy_Ball){
        return heavy_ball();
      }else if(_input.m==Mode::Nesterov){
        return nesterov();
      }
      return gradient_descendant();
    }


    std::vector<value_type> heavy_ball(){

     bool flag=false;
     parameter_type ak=_input.a0;
     unsigned int k=0;
     std::vector<value_type> grad;
     std::vector<value_type> xold(_input.start);
     std::vector<value_type> xnew(_input.start);
     std::vector<value_type> xold1(_input.start);

     if(_input.a==Alpha::Armijo){
      flag=true;
      std::cerr<<"Cannot use Armijo computation for the update of alpha_k and Heavy-Ball method for solving the minimization problem.\n"<<std::endl;
      exit(1);
     }
     while(!flag){
      // update number of iterations
      ++k;

      //gradient evaluation
      grad=compute_grad(xold);

      //xnew computation
      if(k==1){
        xnew=xold-ak*grad;
      }else{
        xnew=xold-ak*grad+_input.eta*(xold-xold1);
      }

      // check convergence
      flag=check_convergence(xold,xnew,k,grad);
      
      // alpha computation
      ak= compute_alpha(xnew,k);

      //update
      xold1=xold;
      xold=xnew;

      }
     std::cout<<"\nNumber of interations: "<<k<<"\n"<<std::endl;
    return xnew;
  }
  // gradient descendant
  std::vector<value_type> gradient_descendant(){

     bool flag=false;
     parameter_type ak=_input.a0;
     unsigned int k=0;
     std::vector<value_type> grad;
     std::vector<value_type> xold(_input.start);
     std::vector<value_type> xnew(_input.start);

     while(!flag){
      // update number of iterations
      ++k;

      //gradient evaluation
      grad=compute_grad(xold);

      //xnew computation
      xnew=xold-ak*grad;

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

// Nesterov method
  std::vector<value_type> nesterov(){

     bool flag=false;
     parameter_type ak=_input.a0;
     unsigned int k=0;
     std::vector<value_type> grad;
     std::vector<value_type>  xold(_input.start);
     std::vector<value_type>  xold1(_input.start);
     std::vector<value_type> xnew(_input.start);

     if(_input.a==Alpha::Armijo){
      flag=true;
      std::cerr<<"Cannot use Armijo computation for the update of alpha_k and Nesterov method for solving the minimization problem.\n"<<std::endl;
      exit(1);
     }
     while(!flag){
      // update number of iterations
      ++k;

      //gradient evaluation
      grad=compute_grad(xold);

      //xnew computation
      if(k==1){
        xnew=xold-ak*grad;
      }else{
        xnew=xold+_input.eta*(xold-xold1)-ak*compute_grad(xold+_input.eta*(xold-xold1));
      }

      // check convergence
      flag=check_convergence(xold,xnew,k,grad);
      
      // alpha computation
      ak=compute_alpha(xnew,k);
      if(1-ak<1){
        set_eta(1-ak);
      }

      //update
      xold1=xold;
      xold=xnew;

      }
     std::cout<<"\nNumber of interations: "<<k<<"\n"<<std::endl;
    return xnew;
  }


    //check convergence
   bool check_convergence(std::vector<value_type> & xold,std::vector<value_type> & xnew, unsigned int k,std::vector<value_type> & grad)const{
        if(k>_input.it){
          std::cout<<"\nLimit of the max number of iterations: no convergence is reached.\n"<<std::endl;
        }else if(euclidean_norm(xold-xnew)<_input.eps_s){
          std::cout<<"\nControl on the step lenght: ?x_k+1 - x_k? < ?s\n"<<std::endl;
        }else if(euclidean_norm(grad)<_input.eps_r){
          std::cout<<"\nControl on the residual: ??f(x_k)? < ?r\n"<<std::endl;
        }
        return k>_input.it or euclidean_norm(xold-xnew)<_input.eps_s or euclidean_norm(grad)<_input.eps_r;
    }


    //compute alpha_k - Armijo rule
   parameter_type Armijo(std::vector<value_type> & xk){
      parameter_type ak=_input.a0;
      value_type gradient_norm_squared=pow(euclidean_norm((evaluate_gradient_diff(xk))),2);
      while(_input.f(xk)-_input.f(xk-ak*evaluate_gradient_diff(xk))<_input.sigma*ak*gradient_norm_squared){
      ak/=2;
    }
    return ak;
    }
    
    //compute alpha_k - exponential decay
    parameter_type exponential_decay(unsigned int k){return _input.a0*exp(-(_input.mu*static_cast<parameter_type>(k)));}

    //compute alpha_k - inverse decay
    parameter_type inverse_decay(unsigned int k){return _input.a0/(1+(_input.mu*static_cast<parameter_type>(k)));}

   
    parameter_type compute_alpha(std::vector<value_type> & xk,unsigned int k){
    if (_input.a == Alpha::Armijo) {
        return Armijo(xk);
    } else if (_input.a == Alpha::exponential_decay) {
        return exponential_decay(k);
    } else if (_input.a == Alpha::inverse_decay) {
        return inverse_decay(k);
    }
    return Armijo(xk);
    }
};

#endif
