#include <iostream>
#include <cmath>
#include <functional>

const double ALPHA = 0.23232;

double f(double x , double t );
void  euler(double &x, double t , double h);

int main (void){
  const double x0 = 1.2345;
  const double h = 0.1;
  const int N = 1000;


  double x = x0 ;
  double t = 0;
  for(int step = 0; step < N; ++step){
    t = 0.0 + step*h;
    std::cout << step*h << "   " << x << std::endl;
    euler(x,t,h);
  }
  
  return  0;
}

double f(double x , double t ){
  return -ALPHA*x;
}

void  euler(double &x, double t , double h){

  x += h*f(x,t);
  
}
