#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"



const double Lx=100, Ly=100;
const int Nx=1,Ny=5,N=Nx*Ny;


const double ZETA=0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double CHI=-0.06626458266981849;

class Cuerpo;
class Colisionador;

class Cuerpo{
private:
  vector3D r,V,F ;double m,R;
public:
  void Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,
	      double m0,double R0);
  void BorreFuerza(void);
  void AgregueFuerza(vector3D F0);
  void Mueva_r(double dt,double Constante);
  void Mueva_V(double dt,double Constante);
  void Dibujese(void);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  friend class Colisionador;
};

void Cuerpo::Inicie(double x0,double y0,double z0,
		    double Vx0,double Vy0,double Vz0,double m0,double R0){
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
  V+=F*(Constante*dt/m); 
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
  std::cout<<"set xrange [-10:110]"<<std::endl;
  std::cout<<"set yrange [-10:110]"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange[0:7]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(void){
  std::cout<<"plot 0,0 ";
  std::cout<<" , "<<100/7<<"*t,0";
  std::cout<<" , "<<100/7<<"*t,100";
  std::cout<<" , 0,"<<100/7<<"*t";
  std::cout<<" , 100,"<<100/7<<"*t";
}
void TermineCuadro(void){
  std::cout<<std::endl;
}
//-------clase colisionador-------

class Colisionador{
private:
  
public:
  void CalculeTodasLasFuerzas(Cuerpo* Grano,double dt);
  void CalculeLaFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2);
			    
};


void Colisionador:: CalculeTodasLasFuerzas(Cuerpo* Grano,double dt){
  //Borrar todas las fuerzas y torques
  for(int ii = 0;ii<N+4;ii++) Grano[ii].BorreFuerza();
  
  //Fuerza entre pares de bolas
  for(int ii =0;ii<N;ii++){
    for(int jj =ii+1;jj<N+4;jj++){
      CalculeLaFuerzaEntre(Grano[ii],Grano[jj]);
    }
  }
}
void  Colisionador::CalculeLaFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2){


}




int main(void){
  double t, dt=1.0e-2;
  double tdibujo;int Ndibujos;
  Cuerpo Grano[N]; 
  Colisionador Newton;
  Crandom ran64(1);
    
  double m0=1, R0=6,V=0,omega0=10,theta;
  double dx=Lx/(Nx+1);double dy=Ly/(Ny+1);
  double T=1, tmax=T;
  
  InicieAnimacion(); Ndibujos=2000;
  //GRANOS
  
  for(int ii=0;ii<Nx;ii++){
    for(int jj=0;jj<Ny;jj++){
      theta=2*M_PI*ran64.r();
      //--------------------(x0       , y0       ,z0, Vx0              , Vy0             ,Vz0, theta0, omega0, m0 , R0);
      Grano[ii+Nx*jj].Inicie((ii+1)*dx, (jj+1)*dy, 0, V*std::cos(theta),V*std::sin(theta),0  , m0 , R0);
    }
  }
  //Newton.CalculeTodasLasFuerzas(Grano,dt);
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    if(tdibujo<tmax/Ndibujos){
      InicieCuadro();
      for(int ii = 0;ii<N;ii++)Grano[ii].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    

    
    //std::cout<<Grano[0].Getx()<<"   "<<Grano[0].Gety()<<std::endl;
    for(int ii = 0;ii<N;ii++)Grano[ii].Mueva_r(dt,ZETA);
    Newton.CalculeTodasLasFuerzas(Grano,dt);
    for(int ii = 0;ii<N;ii++)Grano[ii].Mueva_V(dt,(1-2*LAMBDA)/2);
    for(int ii = 0;ii<N;ii++)Grano[ii].Mueva_r(dt,CHI);
    Newton.CalculeTodasLasFuerzas(Grano,dt);
    for(int ii = 0;ii<N;ii++)Grano[ii].Mueva_V(dt,LAMBDA);
    for(int ii = 0;ii<N;ii++)Grano[ii].Mueva_r(dt,1-2*(CHI+ZETA));
    Newton.CalculeTodasLasFuerzas(Grano,dt);
    for(int ii = 0;ii<N;ii++)Grano[ii].Mueva_V(dt,LAMBDA);		    
    for(int ii = 0;ii<N;ii++)Grano[ii].Mueva_r(dt,CHI);
    Newton.CalculeTodasLasFuerzas(Grano,dt);
    for(int ii = 0;ii<N;ii++)Grano[ii].Mueva_V(dt,(1-2*LAMBDA)/2);
    for(int ii = 0;ii<N;ii++)Grano[ii].Mueva_r(dt,ZETA);
  }
  
  return 0;
}
		    
  
