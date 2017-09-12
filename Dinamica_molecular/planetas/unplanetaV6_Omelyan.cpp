#include <iostream>
#include <cmath>
#include "Vector.h"

const double GM=1.0;
const int N=2;

const double ZETA=0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double CHI=-0.06626458266981849;

class Cuerpo;

class Cuerpo{
private:
  vector3D r,V,F;double m,R;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0, double m0,double R0);
  void CalculeFuerza(void);
  void Mueva_r(double dt,double Constante);
  void Mueva_V(double dt,double Constante);
  void Dibujese(void);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
};

void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0){
  r.cargue(x0,y0,z0);
  V.cargue(Vx0,Vy0,Vz0);
  m=m0;
  R=R0;
}

void Cuerpo::CalculeFuerza(void){
  double aux=GM*m*std::pow(norma2(r),-1.5);
  F=(-aux)*r;
}

void Cuerpo::Mueva_r(double dt,double Constante){
r+=V*(Constante*dt);
}
void Cuerpo::Mueva_V(double dt,double Constante){
V+=F*(Constante*dt);
}

void Cuerpo::Dibujese(void){
  std::cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//------------------Funciones Globales---------

void InicieAnimacion(void){
  //std::cout<<"set terminal gif animate"<<std::endl;
  //std::cout<<"set output 'planeta.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set xrange [-120:120]"<<std::endl;
  std::cout<<"set yrange [-120:120]"<<std::endl;
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
  double t, dt=1.0,tmax;
  double tdibujo;int Ndibujos;
  Cuerpo Planeta;

  double r,omega,V,T,m,R;
  R=5; m=1; r=100; omega=std::sqrt(GM/(r*r*r)); V=omega*r; T=2*M_PI/omega; tmax=1.1*T;  

  
  //-----(x0, y0, Vx0, Vy0, m0, R0);
  //InicieAnimacion(); Ndibujos=50;
  Planeta.Inicie(r, 0, 0, 0, 0.5*V,0, m, R);
  Planeta.CalculeFuerza();
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    /*if(tdibujo<tmax/Ndibujos){
      InicieCuadro();
      Planeta.Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    */
    std::cout<<Planeta.Getx()<<"   "<<Planeta.Gety()<<std::endl;
    Planeta.Mueva_r(dt,ZETA);
    Planeta.CalculeFuerza();
    Planeta.Mueva_V(dt,(1-2*LAMBDA)/2);
    Planeta.Mueva_r(dt,CHI);
    Planeta.CalculeFuerza();
    Planeta.Mueva_V(dt,LAMBDA);
    Planeta.Mueva_r(dt,1-2*(CHI-ZETA));
    Planeta.CalculeFuerza();
    Planeta.Mueva_V(dt,LAMBDA);		    
    Planeta.Mueva_r(dt,CHI);
    Planeta.CalculeFuerza();
    Planeta.Mueva_V(dt,(1-2*LAMBDA)/2);
    Planeta.Mueva_r(dt,ZETA);
  }
  return 0;
}
		    
  
