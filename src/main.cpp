
#include "helper.hpp"
#include "method.hpp"
#include "GetPot"

// function
value_type f1(std::vector<value_type> points){
   return (points[0]*points[1]+4*std::pow(points[0],4)+std::pow(points[1],2)+2*points[0]);
}

// the gradient of the function
std::vector<value_type> df1(std::vector<value_type> points){
   return {points[1]+16*std::pow(points[0],3)+3,points[0]+2*points[1]};
}

int main(int argc,char**argv){
    GetPot command_line(argc,argv);
    input i=read_cl(command_line);
    i.f=f1;
    i.df=df1;
    std::vector<value_type> exact_solution{-0.590551,0.295275};
    i.ex_sol=exact_solution;
    std::vector<value_type> x=solve<Diff::finite_diff,Alpha::Armijo,Mode::gradient>(i);
    std::cout<<"Solution: ";
    print_vector(x);
    return 0;
}