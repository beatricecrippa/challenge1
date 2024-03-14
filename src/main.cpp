#include <functional>
#include <cmath>
#include "method.hpp"
#include "helper.hpp"

value_type f1(vector points){
   return points[0]*points[1]+4*std::pow(points[0],2);
}

vector df1(vector points){
   return {points[1]+6*points[0],points[0]};
}

int main(){
    functionR f1;
    functionRn df1;
    std::vector v={0.0,0.0};
    vector start(v);
    struct input i{f1,df1};
    method m(i);
    vector solution=m.solve();
    return 0;
