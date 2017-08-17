#include<iostream>
#include<cmath>

const double ERR=1.0e-7;
double Omega0=1.0;
double f1(double x1, double x2, double t){
double f1(double x1, double x2, double t){
    
  return x2;
}
double f2(double Omega0,double x1, double x2, double t){
  return -Omega0*Omega0*x1;
}



void UnPasoDeRungeKutta4(double Omega0,double &x1, double &x2, double &t,double dt){
  double dx11,dx21,dx31,dx41;
  double dx12,dx22,dx32,dx42;
  
  dx11=dt*f1(x1,x2,t);                            dx12=dt*f2(Omega0,x1,x2,t);
  dx21=dt*f1(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt);   dx22=dt*f2(Omega0,x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt);
  dx31=dt*f1(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt);   dx32=dt*f2(Omega0,x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt);
  dx41=dt*f1(x1+dx31,x2+0.5*dx32,t+dt);           dx42=dt*f2(Omega0,x1+dx31,x2+dx32,t+dt);
  t+=dt; x1+=(dx11+2*dx21+2*dx31+dx41)/6.0; x2+=(dx12+2*dx22+2*dx32+dx42)/6.0; 
}

double f(double Omega0){
  double x1 = 1.0;
  double x2 = 1.0;
  double t = 0.0;
  double dt=0.01;
  for(int ii = 0; ii<100; ii++){
    // std::cout<<t<<" "<<x1<<" "<<std::endl;
    UnPasoDeRungeKutta4(Omega0,x1,x2,t,dt);
  }
  
  return x1;
  
}

double CerosBiseccion(double a,double b){
  double m;
  double fa,fm;
  //ENTRADA DE DATOS
  fa=f(a);
  //PROCESAMIENTO
  do{ //repita
  // calcule m; calcule fm;
    m=(a+b)/2; fm=f(m);
    if(fa*fm<0)// if(fa y fm tienen signos opuestos)
      {b=m;} //corra b;
  else
    {a=m;fa=fm;}//corra a;
  } while(b-a>ERR);//hasta lograr la precision deseada
  //SALIDA DE DATOS
  return (a+b)/2;
}



int main(void){
  double Omega0 = 0.0;
  for(int ii = 0;ii<10000;ii++){
    std::cout<<Omega0<<" " << f(Omega0) <<std::endl; 
    Omega0+=0.01;
  }
  

  return 0 ;
}
