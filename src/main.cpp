
#include "helper.hpp"
#include "method.hpp"

value_type f1(std::vector<value_type> points){
   return (points[0]*points[1]+4*std::pow(points[0],4)+std::pow(points[1],2)+2*points[0]);
}

std::vector<value_type> df1(std::vector<value_type> points){
   return {points[1]+16*std::pow(points[0],3)+3,points[0]+2*points[1]};
}

int main(int argc, char **argv){
    // argv[1]: Mode (Gradient, Heavy_Ball, Nesterov)
    // argv[2]: Alpha (inverse_decay, exponential_decay, Armijo)
    // argv[3]: Grad computation (User_grad, Finite_diff)
    method m{f1,df1};
    if(argc==4)
{
        m.set_mode(mode_map.at(argv[1]));
        m.set_alpha(alpha_map.at(argv[2]));
        m.set_diffcomp(diff_map.at(argv[3]));
}        
    
   
    std::vector<value_type> solution=m.solve();
    std::cout<<"\nSolution: ";
    print_vector(solution);
    return 0;
}
