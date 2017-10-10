#include <iostream>
#include <cmath>
#include "Vector.h"

const double epsilon = 1.0;
const double r0 = 10;
const double KbT = 0.5;
const int N=1;

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
  double GetVx(void){return V.x();};
  double GetV(void){return norma(V);};
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
V+=F*(Constante*dt)/m;
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
  std::cout<<"set xrange [-0:20]"<<std::endl;
  std::cout<<"set yrange [-10:10]"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange[0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(void){
  std::cout<<"plot 0,0 ";
  //std::cout<<" , "<<100/7<<"*t,0";
  //std::cout<<" , "<<100/7<<"*t,100";
  //std::cout<<" , 0,"<<100/7<<"*t";
  //std::cout<<" , 100,"<<100/7<<"*t";
}
void TermineCuadro(void){
  std::cout<<std::endl;
}
//-------clase colisionador-------

class Colisionador{
private:

public:
 void  CalculeTodasLasFuerzas(Cuerpo* Particula);
 void  CalculeLaFuerzaEntre(Cuerpo & Particula1,Cuerpo & Particula2);
  
};

void Colisionador:: CalculeTodasLasFuerzas(Cuerpo* Particula){
  //Borrar todas las fuerzas
  for(int ii = 0;ii<N;ii++) Particula[ii].BorreFuerza();
  for(int ii =0;ii<N;ii++){
    CalculeLaFuerzaEntre(Particula[ii],Particula[N]);
  }
}
void  Colisionador::CalculeLaFuerzaEntre(Cuerpo & Particula1,Cuerpo & Particula2){
  vector3D F2, dr, r_unitario; double Normadr,F_LJ;
  dr=Particula1.r-Particula2.r ; Normadr=norma(dr);r_unitario=dr/Normadr;
  //calculo fuerza de Lennard-Jones 
  F_LJ=12*epsilon/Normadr*(std::pow(r0/Normadr,12)-std::pow(r0/Normadr,6));
  F2=F_LJ*r_unitario;
  Particula1.AgregueFuerza(F2);
  
}



int main(void){
  double t, dt=1.0e-2;
  double tdibujo;int Ndibujos;
  Cuerpo Particula[N+1]; 
  Colisionador Newton;
  
    
  double m0=1.0, R0=2.5,Vx0=std::sqrt(2*KbT/m0),Vy0=0,x0=10,y0=0;
  double  tmax=100;

  //InicieAnimacion();
  Ndibujos=2000;

  //PARTICULA
  //-------------(x0, y0,z0, Vx0, Vy0,Vz0, m0, R0);
  for(int ii=0;ii<N;ii++){
    Particula[ii].Inicie(x0, y0, 0, Vx0, Vy0,  0, m0, R0);
  }
  //PARTICULA CENTRAL
  Particula[N].Inicie(0, 0, 0, 0, 0,  0, 1, 1);
  
  Newton.CalculeTodasLasFuerzas(Particula);
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    /*if(tdibujo<tmax/Ndibujos){
      InicieCuadro();
      for(int ii = 0;ii<N;ii++)Particula[ii].Dibujese();
      TermineCuadro();
      tdibujo=0;
      }
    */
    
    
    std::cout<<t<<"   "<<Particula[0].Getx()<<std::endl;
    for(int ii = 0;ii<N;ii++)Particula[ii].Mueva_r(dt,ZETA);
    Newton.CalculeTodasLasFuerzas(Particula);
    for(int ii = 0;ii<N;ii++)Particula[ii].Mueva_V(dt,(1-2*LAMBDA)/2);
    for(int ii = 0;ii<N;ii++)Particula[ii].Mueva_r(dt,CHI);
    Newton.CalculeTodasLasFuerzas(Particula);
    for(int ii = 0;ii<N;ii++)Particula[ii].Mueva_V(dt,LAMBDA);
    for(int ii = 0;ii<N;ii++)Particula[ii].Mueva_r(dt,1-2*(CHI+ZETA));
    Newton.CalculeTodasLasFuerzas(Particula);
    for(int ii = 0;ii<N;ii++)Particula[ii].Mueva_V(dt,LAMBDA);		    
    for(int ii = 0;ii<N;ii++)Particula[ii].Mueva_r(dt,CHI);
    Newton.CalculeTodasLasFuerzas(Particula);
    for(int ii = 0;ii<N;ii++)Particula[ii].Mueva_V(dt,(1-2*LAMBDA)/2);
    for(int ii = 0;ii<N;ii++)Particula[ii].Mueva_r(dt,ZETA);
  }
  
return 0;
}
		    
  
