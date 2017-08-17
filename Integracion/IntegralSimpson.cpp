#include<iostream>
#include<fstream> //Para leer archivos
#include <cmath>
using namespace std; 

const double ERR=1e-7;

double f(double alpha, double x, double t){
  return cos(alpha*t-x*sin(t));
}

double IntegralPorSimpson(double alpha, double x, double a, double b, int N){ //N corresponde a la particion entre a y b, ademas debe ser par
  double h = (b-a)/N, suma, t;
  int i;
  for(suma=0, i=0; i<=N; i++ ){
    double t= a+h*i;
    if( i==0 || i==N)
      suma+=f(alpha,x,t);
    else if (i%2==0)
      suma+=2*f(alpha,x,t);
    else
      suma+=4*f(aĺpha,x,t);
  }

  return suma*h/3.0;
}
double Bessel(double alpha, double x){
  return IntegralPorSimpson(double alpha, double x, 0, M_PI, 50)/M_PI;
}

int main(void){
  double x,alpha=0;
  for(x=0;N<20;N+=0.1){
    h=(b-a)/N;
    cout<< x << " " << Bessel(alpha,x)-2<<endl; 
  }
  return 0;
}
