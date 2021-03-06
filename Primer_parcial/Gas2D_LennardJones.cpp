#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"

const double epsilon = 1.0;
const double r0 = 10;
const double KbT =0.1;

const double K=50;
const double Lx=100, Ly=100;
const int Nx=5,Ny=5,N=Nx*Ny;

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
  std::cout<<"set xrange [-40:100]"<<std::endl;
  std::cout<<"set yrange [-40:100]"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange[0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(void){
  std::cout<<"plot 0,0 ";
  // std::cout<<" , "<<100/7<<"*t,0";
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
  void  CalculeLaFuerzaEntreParticulas(Cuerpo & Particula1,Cuerpo & Particula2);
  void  CalculeLaFuerzaEntreParticulasYParedes(Cuerpo & Particula1,Cuerpo & Particula2);
};

void Colisionador:: CalculeTodasLasFuerzas(Cuerpo* Particula){
  //Borrar todas las fuerzas
  for(int ii = 0;ii<N+4;ii++) Particula[ii].BorreFuerza();
  for(int ii =0;ii<N;ii++){
    for(int jj =ii+1;jj<N;jj++){
      CalculeLaFuerzaEntreParticulas(Particula[ii],Particula[jj]);
    }
  }

  for(int ii =0;ii<N;ii++){
    CalculeLaFuerzaEntreParticulasYParedes(Particula[ii],Particula[N]);
    CalculeLaFuerzaEntreParticulasYParedes(Particula[ii],Particula[N+1]);
    CalculeLaFuerzaEntreParticulasYParedes(Particula[ii],Particula[N+2]);
    CalculeLaFuerzaEntreParticulasYParedes(Particula[ii],Particula[N+3]);
  }
}
void  Colisionador::CalculeLaFuerzaEntreParticulas(Cuerpo & Particula1,Cuerpo & Particula2){
  vector3D F2, dr, r_unitario; double Normadr,F_LJ;
  dr=Particula1.r-Particula2.r ; Normadr=norma(dr);r_unitario=dr/Normadr;
  //calculo fuerza de Lennard-Jones 
  F_LJ=12*epsilon/Normadr*(std::pow(r0/Normadr,12)-std::pow(r0/Normadr,6));
  F2=F_LJ*r_unitario;
  Particula1.AgregueFuerza(F2);Particula2.AgregueFuerza((-1)*F2);
}
void  Colisionador::CalculeLaFuerzaEntreParticulasYParedes(Cuerpo & Particula1,Cuerpo & Particula2){
  vector3D F2, dr, r_unitario; double Normadr,s;
  dr=Particula2.r-Particula1.r; Normadr=norma(dr);r_unitario=dr/Normadr;
  s=(Particula1.R+Particula2.R)-Normadr;
  
  if(s>0){//si se chocan
    F2=r_unitario*(K*std::pow(s,1.5));
    Particula1.AgregueFuerza(F2*(-1));   Particula2.AgregueFuerza(F2);
  }
  
}



int main(void){
  double t, dt=1.0e-3;
  double tdibujo;int Ndibujos;
  Cuerpo Particula[N+4]; 
  Colisionador Newton;
  Crandom ran64(1);double theta;
    
  double m0=1.0, R0=2.5,V=std::sqrt(2*KbT/m0);
  double Rpared=10000,Mpared=1000*m0;
  double dx=Lx/(Nx+1);double dy=Ly/(Ny+1);
  double tmax=10;

  InicieAnimacion();
  Ndibujos=50;
  //PAREDES
  //---------------(x0       , y0      ,z0 , Vx0, Vy0,Vz0 , m0    , R0);
  //Pared arriba
  Particula[N  ].Inicie(70    ,100+Rpared,  0, 0  , 0,0, Mpared, Rpared);
  //Pared abajo
  Particula[N+1].Inicie(70     ,-Rpared-40  ,  0, 0  , 0,0, Mpared, Rpared);
  //Pared derecha
  Particula[N+2].Inicie(100+Rpared, 70    ,  0, 0  , 0,0, Mpared, Rpared);
  //Pared izquerda
  Particula[N+3].Inicie(-Rpared-40  , 70    ,  0, 0  , 0,0, Mpared, Rpared);

  //GRANOS
  //-------------(x0  , y0  ,z0, Vx0  , Vy0 , m0    , R0);
  for(int ii=0;ii<Nx;ii++){
    for(int jj=0;jj<Ny;jj++){
      theta=2*M_PI*ran64.r();
      Particula[ii+Nx*jj].Inicie((ii+1)*10, (jj+1)*10, 0, V*std::cos(theta),V*std::sin(theta),0, m0, R0);
    }
  }
  Newton.CalculeTodasLasFuerzas(Particula);
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    if(tdibujo<tmax/Ndibujos){
      InicieCuadro();
      for(int ii = 0;ii<N;ii++)Particula[ii].Dibujese();
      TermineCuadro();
      tdibujo=0;
     }
    

    
    //std::cout<<Particula[0].Getx()<<"   "<<Particula[0].Gety()<<std::endl;
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
  //for(int ii = 0;ii<N;ii++)std::cout<<Particula[ii].GetVx()<<std::endl;
  return 0;
}
		    
  
