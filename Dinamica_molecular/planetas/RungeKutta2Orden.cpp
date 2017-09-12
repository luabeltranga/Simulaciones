#include<iostream>
#include<fstream> 
#include <cmath>

const double ERR=1e-7;

double f1(double x1, double x2, double t);
double f2(double Omega0, double x1, double x2, double t);
void UnPasoDeRungeKutta4(double Omega0, double & x1, double & x2, double & t, double & dt);
double f(double Omega0);
double CerosBiseccionDefEntre(double a,double b);


int main(void){


  double a=2,b=5;
  double Omega = 0.0;
  for(int ii = 0; ii<100;ii++){
    std::cout<<Omega<<" "<<f(Omega)<<std::endl;
    Omega+=0.1;
  }
  //cout <<"El cero entre a="<<a<<" y b="<<b<<" es "<<CerosBiseccionDefEntre(a,b)<<endl;
  return 0;


}

double f1(double x1, double x2, double t){
  return x2;
}

double f2(double Omega0, double x1, double x2, double t){
  return -Omega0*Omega0*x1;
}

void UnPasoDeRungeKutta4(double Omega0, double & x1, double & x2, double & t, double & dt){ 
  double dx11,dx21,dx31,dx41;                   double dx12,dx22,dx32,dx42;
  dx11=dt*f1(x1,x2,t);                          dx12=dt*f2(Omega0,x1,x2,t);
  dx21=dt*f1(x1+dx11/2,x2+dx12/2,t+dt/2);       dx22=dt*f2(Omega0,x1+dx11/2,x2+dx12/2,t+dt/2);
  dx31=dt*f1(x1+dx21/2,x2+dx22/2,t+dt/2);       dx32=dt*f2(Omega0,x1+dx21/2,x2+dx22/2,t+dt/2);
  dx41=dt*f1(x1+dx31,x2+dx32,t+dt);             dx42=dt*f2(Omega0,x1+dx31,x2+dx32,t+dt);                          
  t+=dt; x1+=(dx11+2*dx21+2*dx31+dx41)/6;       x2+=(dx12+2*dx22+2*dx32+dx42)/6;
}

double f(double Omega0){
  double x1=0,x2=1,t=0;
  double dt=0.01;
  for(t=0;t<1; ){
    UnPasoDeRungeKutta4(Omega0,x1,x2,t,dt);
  }
  return x1;
}

double CerosBiseccionDefEntre(double a,double b){
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



