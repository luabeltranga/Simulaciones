#include <iostream>
#include <cmath>
#include "Random64.h"

const double Deltat=0.01; //en picosegunods
const double L=100; //en Amstrong
const double T=300; //en kelvin
const double MASA=22.8916; //en u.m.a
const double e=1; //carga electron
const double D=0.132; //en A²/ps
const double kb=0.826; //en esas unidades
const double Gamma=kb*T/(MASA*D); 
const double sigma=std::sqrt(2*D*Deltat);
const double dtUmGamma=Deltat/(2*MASA*Gamma);
const int N=4000;

class Cuerpo;

class Cuerpo{
private:
  double x,Fx,Fxold,m,R;
public:
  void Inicie(double x0,double m0,double R0);
  void Muevase(Crandom &ran64);
  void CalculeFuerza(double E);
  void Dibujese(void);
  double Getx(void){return x;};
};

void Cuerpo::Inicie(double x0,double m0,double R0){
  x=x0;
  m=m0;
  R=R0;
  Fx=0;
}

void Cuerpo::CalculeFuerza(double E){
  Fxold=Fx;
  Fx=e*E;

}

void Cuerpo::Muevase(Crandom & ran64){
  x+=dtUmGamma*(3*Fx-Fxold)+ran64.gauss(0,sigma);

}
void Cuerpo::Dibujese(void){
  std::cout<<", "<<x<<"+"<<R<<"*cos(t),"<<0<<"+"<<R<<"*sin(t)";
}
//------------------Funciones Globales---------

void InicieAnimacion(void){
  //std::cout<<"set terminal gif animate"<<std::endl;
  //std::cout<<"set output 'planeta.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set xrange [-10:110]"<<std::endl;
  std::cout<<"set yrange [-10:10]"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange[0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(void){
  std::cout<<"plot 0,0 ";
}
void TermineCuadro(void){
  std::cout<<std::endl;
}


int main(void){
  double t; double tmax=40;
  Cuerpo Na[N];
  Crandom ran64(1);
  double E=0;
  double sigmax0=5;
  
  InicieAnimacion(); 
  //----------( x0,   m0, R0);
  for(int ii=0;ii<N;++ii) Na[ii].Inicie(ran64.gauss(L/2,sigmax0), MASA, 4);
  for(int ii=0;ii<N;++ii) Na[ii].CalculeFuerza(E);
  for(t=0;t<tmax;t+=Deltat){
    InicieCuadro();
    Na[0].Dibujese();
    TermineCuadro();
    
    for(int ii=0;ii<N;++ii) Na[ii].CalculeFuerza(E);
    for(int ii=0;ii<N;++ii) Na[ii].Muevase(ran64);
  }

  //for(int ii=0;ii<N;++ii) std::cout<<Na[ii].Getx()<<std::endl;
  
  return 0;
}




