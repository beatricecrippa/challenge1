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
    method()=default;
    method(input &i): _input(i){print_struct(_input);};
    //constructor by passing as parameter only f and df
    method(functionR  _f, functionRn  _df){
        _input.f=_f;
        _input.df=_df;
    };
     void set_mode(Mode _m){
         _input.m=_m;
     }

     void set_alpha(Alpha _a){
         _input.a=_a;
     }
     void set_diffcomp(Diff _d){
         _input.d=_d;
     }
    void set_eta(parameter_type value){
      _input.eta=value;
    }

    
    // solving the minimization problem
    std::vector<value_type> solve(){
      print_struct(_input);
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
     std::vector<value_type> grad((_input.start).size(),0);
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
      switch(_input.d){
        case Diff::Finite_diff:
         for(size_t i=0;i<xold.size();++i){
          std::vector<value_type> v1(xold),v2(xold);
          v1[i]+=h;
          v2[i]-=h;
          grad[i]=((_input.f(v1)-_input.f(v2))/(2.*h));
         }
         break;
         case Diff::User_grad:
         grad=_input.df(xold);
         break;
      }

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
     std::cout<<"Number of interations: "<<k<<"\n"<<std::endl;
    return xnew;
  }
  // gradient descendant
  std::vector<value_type> gradient_descendant(){

     bool flag=false;
     parameter_type ak=_input.a0;
     unsigned int k=0;
      std::vector<value_type> grad((_input.start).size(),0);
     std::vector<value_type> xold(_input.start);
     std::vector<value_type> xnew(_input.start);

     while(!flag){
      // update number of iterations
      ++k;

      //gradient evaluation
      switch(_input.d){
        case Diff::Finite_diff:
         for(size_t i=0;i<xold.size();++i){
          std::vector<value_type> v1(xold),v2(xold);
          v1[i]+=h;
          v2[i]-=h;
          grad[i]=((_input.f(v1)-_input.f(v2))/(2.*h));
         }
         break;
         case Diff::User_grad:
         grad=_input.df(xold);
         break;
      }

      //xnew computation
      xnew=xold-ak*grad;

      // check convergence
      flag=check_convergence(xold,xnew,k,grad);
      
      // alpha computation
      ak= compute_alpha(xnew,k);

      //update
      xold=xnew;

      }
     std::cout<<"Number of interations: "<<k<<"\n"<<std::endl;
    return xnew;
  }

// Nesterov method
  std::vector<value_type> nesterov(){

     bool flag=false;
     parameter_type ak=_input.a0;
     unsigned int k=0;
      std::vector<value_type> grad((_input.start).size(),0);
      std::vector<value_type> grad_y((_input.start).size(),0);
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
      switch(_input.d){
        case Diff::Finite_diff:
         for(size_t i=0;i<xold.size();++i){
          std::vector<value_type> v1(xold),v2(xold);
          v1[i]+=h;
          v2[i]-=h;
          grad.emplace_back((_input.f(v1)-_input.f(v2))/(2.*h));
         }
         break;
         case Diff::User_grad:
         grad=_input.df(xold);
         break;
      }

        //gradient y evaluation
      switch(_input.d){
        case Diff::Finite_diff:
         for(size_t i=0;i<(xold+_input.eta*(xold-xold1)).size();++i){
          std::vector<value_type> v1(xold+_input.eta*(xold-xold1)),v2(xold+_input.eta*(xold-xold1));
          v1[i]+=h;
          v2[i]-=h;
          grad_y[i]=(_input.f(v1)-_input.f(v2))/(2.*h);
         }
         break;
         case Diff::User_grad:
         grad_y=_input.df(xold+_input.eta*(xold-xold1));
         break;
      }

      //xnew computation
      if(k==1){
        xnew=xold-ak*grad;
      }else{
        xnew=xold+_input.eta*(xold-xold1)-ak*grad_y;
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
     std::cout<<"Number of interations: "<<k<<"\n"<<std::endl;
    return xnew;
  }


    //check convergence
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


    //compute alpha_k - Armijo rule
   parameter_type Armijo(std::vector<value_type> & xk){
      parameter_type ak=_input.a0;
       std::vector<value_type> grad((_input.start).size(),0);
      switch(_input.d){
        case Diff::Finite_diff:
         for(size_t i=0;i<xk.size();++i){
          std::vector<value_type> v1(xk),v2(xk);
          v1[i]+=h;
          v2[i]-=h;
          grad[i]=(_input.f(v1)-_input.f(v2))/(2.*h);
         }
         break;
         case Diff::User_grad:
         grad=_input.df(xk);
         break;
      }

      value_type gradient_norm_squared=pow(euclidean_norm((grad)),2);
      while(_input.f(xk)-_input.f(xk-ak*grad)<_input.sigma*ak*gradient_norm_squared){
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
