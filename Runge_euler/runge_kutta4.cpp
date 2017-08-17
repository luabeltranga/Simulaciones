#include<iostream>
#include<cmath>

double f(double x,double t){
  return x;
}
double fexact(double t){
  return std::exp(t);
}

void UnPasoDeRungeKutta4(double &x, double &t,double dt){
  double dx1,dx2,dx3,dx4;
  dx1=dt*f(x,t);
  dx2=dt*f(x+0.5*dx1,t+0.5*dt);
  dx3=dt*f(x+0.5*dx2,t+0.5*dt);
  dx4=dt*f(x+dx3,t+dt);
  t+=dt; x+=(dx1+2*dx2+2*dx3+dx4)/6.0; 
}

int main(void){
  double x = 1.0;
  double t = 0.0;
  double dt=0.01;
  for(int ii = 0; ii<1000; ii++){
    std::cout<<t<<" "<<x<<" "<<fexact(t)<<std::endl;
    UnPasoDeRungeKutta4(x,t,dt);
  }

  return 0 ;
}
