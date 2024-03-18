
#include "helper.hpp"
#include "method.hpp"
#include "GetPot"

value_type f1(vector points){
   return (points[0]*points[1]+4*std::pow(points[0],4)+std::pow(points[1],2)+2*points[0]);
}

vector df1(vector points){
   return {points[1]+16*std::pow(points[0],3)+3,points[0]+2*points[1]};
}

int main(int argc,char ** argv){
   GetPot command_line(argc,argv);
   
    method m{f1,df1};
    vector solution=m.solve();
    print(solution);
    return 0;
}
