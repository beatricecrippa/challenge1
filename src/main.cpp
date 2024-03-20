
#include "helper.hpp"
#include "method.hpp"

value_type f1(std::vector<value_type> points){
   return (points[0]*points[1]+4*std::pow(points[0],4)+std::pow(points[1],2)+2*points[0]);
}

std::vector<value_type> df1(std::vector<value_type> points){
   return {points[1]+16*std::pow(points[0],3)+3,points[0]+2*points[1]};
}

int main(int argc, char **argv){
    // argv[1]: Mode (gradient, Heavy_Ball, Nesterov)
    // argv[2]: Alpha (inverse_decay, exponential_decay, Armijo)
    // argv[3]: Grad computation (user_grad, finite_diff)
    input i;
    i.f=f1;
    i.df=df1;
    read_inputs(argv,argc,i);
    method m(i);
    std::vector<value_type> solution=m.solve();
    std::cout<<"Solution: ";
    print_vector(solution);
    return 0;
}
