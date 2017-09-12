#include <iostream>
#include <cmath>
#include "Vector.h"

const double G=1.0;
const int N=3;
const double dt=0.1;

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
  void Dibujese_Rotado(double t);
  void Datos_Rotado(double t, double omega);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  double GetVx(void){return V.x();};
  double GetVy(void){return V.y();};
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
void Cuerpo::Dibujese_Rotado(double t){
  double theta = norma(V)*t/norma(r);
  double rx_rotado = r.x()*std::cos(theta)+r.y()*std::sin(theta);
  double ry_rotado = -r.y()*std::cos(theta)+r.x()*std::sin(theta);
  std::cout<<", "<<rx_rotado<<"+"<<R<<"*cos(t),"<<ry_rotado<<"+"<<R<<"*sin(t)";
}
void Cuerpo::Datos_Rotado(double t, double omega){
  double theta=omega*t;
  double rx_rotado = r.x()*std::cos(theta)+r.y()*std::sin(theta);
  double ry_rotado = -r.y()*std::cos(theta)+r.x()*std::sin(theta);
  std::cout<<t<<"    "<<rx_rotado<<"    ";
}
//------------------Funciones Globales---------

void InicieAnimacion(void){
  //std::cout<<"set terminal pdf enhanced"<<std::endl;
  //std::cout<<"set output 'planeta.pdf'"<<std::endl;
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
    for(int jj =ii+1;jj<N;jj++){
      CalculeLaFuerzaEntre(Planeta[ii],Planeta[jj]);
    }
  }
}
void  Colisionador::CalculeLaFuerzaEntre(Cuerpo & Planeta1,Cuerpo & Planeta2){
  vector3D dr=Planeta2.r-Planeta1.r;
  vector3D F1;
  double aux=G*Planeta1.m*Planeta2.m*std::pow(norma2(dr),-1.5);
  F1=dr*aux;
  Planeta1.AgregueFuerza(F1);   Planeta2.AgregueFuerza(-1*F1);
}



int main(void){
  
  double t;
  double tdibujo;int Ndibujos;
  Cuerpo Planeta[N]; 
  Colisionador Newton;
  //masa sol
  double m0=1047;
  //masa jupiter
  double m1=1;
  //masa troyano
  double m2 = 0.005;
  //distancia entre el sol y jupiter
  double r=1000;

  double R0=100, R1=50, R2=50;

  double M=m0+m1+m2;

  double x0=-(m1/M)*r, x1=r+x0;

  double omega=std::sqrt(G*M/(r*r*r)), Vy0=omega*x0, Vy1=omega*x1, T=2*M_PI/omega, tmax=20*T;  

  //InicieAnimacion(); Ndibujos=500;
  //-----(x0, y0, Vx0, Vy0, m0, R0);
  Planeta[0].Inicie(x0, 0, 0, 0, Vy0,0, m0, R0);
  Planeta[1].Inicie(x1, 0, 0, 0, Vy1,0, m1, R1);
  Planeta[2].Inicie(0.5*x1, -std::sqrt(3)*x1/2, 0, std::sqrt(3)*Vy1/2, 0.5*Vy1,0, m2, R2);

  Newton.CalculeTodasLasFuerzas(Planeta);
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    if(tdibujo<tmax/Ndibujos){
      //InicieCuadro();
      //for(int ii = 0;ii<N;ii++)Planeta[ii].Dibujese_Rotado(t);
      for(int ii = 0;ii<N;ii++)Planeta[ii].Datos_Rotado(t,omega);
      TermineCuadro();
      tdibujo=0;
    }
    
    //double r=std::pow(Planeta[1].Getx(),2)+std::pow(Planeta[1].Gety(),2);
    //double V=std::pow(Planeta[1].GetVx(),2)+std::pow(Planeta[1].GetVy(),2);
    //std::cout<<omega*omega<<"  "<<V/r<<std::endl;
    
    //std::cout<<Planeta[0].Getx()<<"   "<<Planeta[0].Gety()<<std::endl;
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
  }
  
  return 0;
}
		    
  
