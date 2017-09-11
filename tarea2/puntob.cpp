
#include <iostream>
#include <cmath>
#include "Vector.h"
#include <fstream>
const double G=1.0;
const int N=2;

const double ZETA=0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double CHI=-0.06626458266981849;

class Cuerpo;
class Colisionador;

class Cuerpo{
private:
  vector3D r,V,F;double m,R;
public:
  void Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0, double m0,double R0);
  void BorreFuerza(void);
  void AgregueFuerza(vector3D F0);
  void Mueva_r(double dt,double Constante);
  void Mueva_V(double dt,double Constante);
  void Dibujese(void);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  friend class Colisionador;
};

void Cuerpo::Inicie(double x0,double y0,double z0,double Vx0,double Vy0,double Vz0,double m0,double R0){
  r.cargue(x0,y0,z0);
  V.cargue(Vx0,Vy0,Vz0);
  m=m0;
  R=R0;
}

void Cuerpo::BorreFuerza(void){
  F.cargue(0,0,0);
}

void Cuerpo::AgregueFuerza(vector3D F0){
  F+=F0;
}

void Cuerpo::Mueva_r(double dt,double Constante){
r+=V*(Constante*dt);
}
void Cuerpo::Mueva_V(double dt,double Constante){
V+=F/m*(Constante*dt);
}

void Cuerpo::Dibujese(void){
  std::cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//------------------Funciones Globales---------

void InicieAnimacion(void){
  //std::cout<<"set terminal gif animate delay 10"<<std::endl;
  //std::cout<<"set output 'planeta.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set xrange [-2000:2000]"<<std::endl;
  std::cout<<"set yrange [-2000:2000]"<<std::endl;
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
//-------clase colisionador-------

class Colisionador{
private:

public:
 void  CalculeTodasLasFuerzas(Cuerpo* Planeta);
 void  CalculeLaFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2);
  
};

void Colisionador:: CalculeTodasLasFuerzas(Cuerpo* Planeta){
  //Borrar todas las fuerzas
  for(int ii = 0;ii<N;ii++) Planeta[ii].BorreFuerza();
  for(int ii =0;ii<N;ii++){
    for(int jj =0;jj<ii;jj++){
      CalculeLaFuerzaEntre(Planeta[ii],Planeta[jj]);
    }
  }
}
void  Colisionador::CalculeLaFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2){
  vector3D dr=Planeta2.r-Planeta1.r;
  vector3D F2;
  double aux=-G*Planeta1.m*Planeta2.m*std::pow(norma2(dr),-1.5);
  F2=dr*aux;
  Planeta2.AgregueFuerza(F2);   Planeta1.AgregueFuerza(F2*(-1));
}



int main(void){
  
  double t, dt=0.1;
  double tdibujo;int Ndibujos;
  Cuerpo Planeta[N]; 
  Colisionador Newton;
  
  double m0=1047,m1=1,r=1000;

  double R0=30, R1=10;

  double M=m0+m1;

  double x0=-(m1/M)*r, x1=x0+r;

  double omega=std::sqrt(G*M/(r*r*r)), Vy0=omega*x0, Vy1=omega*x1, T=2*M_PI/omega, tmax=20*T;  
  //std::cout<<omega<<std::endl;
  ofstream jupiter;
  
  //-----(x0, y0, Vx0, Vy0, m0, R0);
  //sol
  Planeta[0].Inicie(x0, 0, 0, 0, Vy0,0, m0, R0);
  //jupiter
  Planeta[1].Inicie(x1, 0, 0, 0, Vy1,0, m1, R1);
  jupiter.open("jupiter.dat");
  Newton.CalculeTodasLasFuerzas(Planeta);
  
  jupiter<<Planeta[1].Getx()<<"   "<<Planeta[1].Gety()<<std::endl;
  for(t=dt;t<tmax;t+=dt){
    for(int ii = 0;ii<N;ii++)Planeta[ii].Mueva_r(dt,ZETA);
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(int ii = 0;ii<N;ii++)Planeta[ii].Mueva_V(dt,(1-2*LAMBDA)/2);
    for(int ii = 0;ii<N;ii++)Planeta[ii].Mueva_r(dt,CHI);
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(int ii = 0;ii<N;ii++)Planeta[ii].Mueva_V(dt,LAMBDA);
    for(int ii = 0;ii<N;ii++)Planeta[ii].Mueva_r(dt,1-2*(CHI+ZETA));
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(int ii = 0;ii<N;ii++)Planeta[ii].Mueva_V(dt,LAMBDA);		    
    for(int ii = 0;ii<N;ii++)Planeta[ii].Mueva_r(dt,CHI);
    Newton.CalculeTodasLasFuerzas(Planeta);
    for(int ii = 0;ii<N;ii++)Planeta[ii].Mueva_V(dt,(1-2*LAMBDA)/2);
    for(int ii = 0;ii<N;ii++)Planeta[ii].Mueva_r(dt,ZETA);
    jupiter<<Planeta[1].Getx()*std::cos(omega*t)+Planeta[1].Gety()*std::sin(omega*t);
    jupiter<<"   "<<Planeta[1].Gety()*std::cos(omega*t)-Planeta[1].Getx()*std::sin(omega*t);
    jupiter<<"   "<<Planeta[0].Getx()*std::cos(omega*t)+Planeta[0].Gety()*std::sin(omega*t);
    jupiter<<"   "<<Planeta[0].Gety()*std::cos(omega*t)-Planeta[0].Getx()*std::sin(omega*t)<<std::endl;
  }
  jupiter.close();
    
  return 0;
}
		    
  
