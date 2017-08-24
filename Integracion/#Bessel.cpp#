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
    t= a+h*i;
    if( i==0 || i==N)
      suma+=f(alpha,x,t);
    else if (i%2==0)
      suma+=2*f(alpha,x,t);
    else
      suma+=4*f(alpha,x,t);
  }
 
  return suma*h/3.0;
}
double Bessel(double alpha, double x){
  return IntegralPorSimpson(alpha,x,0,M_PI,50)/M_PI;
}

double CerosBiseccion(double a,double b,double alpha){
  double m;
  double fa,fm;
  //ENTRADA DE DATOS
  fa=Bessel(alpha,a);
  //PROCESAMIENTO
  do{ //repita
  // calcule m; calcule fm;
    m=(a+b)/2; fm=Bessel(alpha,m);
    if(fa*fm<0)// if(fa y fm tienen signos opuestos)
      {b=m;} //corra b;
  else
    {a=m;fa=fm;}//corra a;
  } while(b-a>ERR);//hasta lograr la precision deseada
  //SALIDA DE DATOS
  return (a+b)/2;
}

int main(void){
  double x,alpha=0;
  //for(x=0;x<20;x+=0.1){
  //cout<< x << " " << Bessel(alpha,x)<<endl;
  //}
  double a=0;
  double b=1;
  double c,fa,fb;
  for(int i=0; i<20;i++){
    a+=1;
    b+=1;
    fa=Bessel(alpha,a),fb=Bessel(alpha,b);
    c=CerosBiseccion(a,b,alpha);
    //cout << c << fa <<fb<< endl;
    if(fa*fb<0)
      {cout<<"El cero esta en x="<< c <<endl;}
    // else
    //{cout<<"No hay cero"<<endl;}
  }
  return 0;
}
