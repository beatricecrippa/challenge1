
#include "helper.hpp"
#include "method.hpp"

// function
value_type f1(std::vector<value_type> points){
   return (points[0]*points[1]+4*std::pow(points[0],4)+std::pow(points[1],2)+2*points[0]);
}

// the gradient of the function
std::vector<value_type> df1(std::vector<value_type> points){
   return {points[1]+16*std::pow(points[0],3)+3,points[0]+2*points[1]};
}

int main(){
    input i;
    i.f=f1;
    i.df=df1;

    std::vector<value_type> solution=solve<Diff::finite_diff,Alpha::inverse_decay,Mode::gradient>(i);
    std::cout<<"Solution: ";
    print_vector(solution);
    return 0;
}
