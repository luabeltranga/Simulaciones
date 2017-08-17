#include<iostream>
#include<cmath>

const double Omega0=1.0;

double f1(double x1, double x2, double t){
  return x2;
}
double f2(double x1, double x2, double t){
  return -Omega0*Omega0*x1;
}
double fexact(double t){
  return std::exp(t);
}

void UnPasoDeRungeKutta4(double &x1, double &x2, double &t,double dt){
  double dx11,dx21,dx31,dx41;
  double dx12,dx22,dx32,dx42;
  
  dx11=dt*f1(x1,x2,t);                            dx12=dt*f2(x1,x2,t);
  dx21=dt*f1(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt);   dx22=dt*f2(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt);
  dx31=dt*f1(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt);   dx32=dt*f2(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt);
  dx41=dt*f1(x1+dx31,x2+0.5*dx32,t+dt);           dx42=dt*f2(x1+dx31,x2+dx32,t+dt);
  t+=dt; x1+=(dx11+2*dx21+2*dx31+dx41)/6.0; x2+=(dx12+2*dx22+2*dx32+dx42)/6.0; 
}

int main(void){
  double x1 = 1.0;
  double x2 = 1.0;
  double t = 0.0;
  double dt=0.01;
  for(int ii = 0; ii<1000; ii++){
    std::cout<<t<<" "<<x1<<" "<<std::endl;
    UnPasoDeRungeKutta4(x1,x2,t,dt);
  }

  return 0 ;
}
