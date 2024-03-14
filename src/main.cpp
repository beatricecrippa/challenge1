#include <functional>
#include <cmath>
#include "helper.hpp"
#include "method.hpp"

value_type f1(vector points){
   return points[0]*points[1]+4*std::pow(points[0],4)+std::pow(points[1],2)+2*points[0];
}

vector df1(vector points){
   return {points[1]+16*std::pow(points[0],3)+3,points[0]+2*points[1]};
}

int main(){
    vector start({0.0,0.0});
    method m{f1,df1};
    vector solution=m.solve();
    return 0;
}
