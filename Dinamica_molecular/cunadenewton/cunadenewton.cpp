#include <iostream>
#include <cmath>
#include "Vector.h"

const double K=1.0e4;
const double g=9.8;
const int N=3;

const double ZETA=0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double CHI=-0.06626458266981849;

class Cuerpo;
class Colisionador;

class Cuerpo{
private:
  double tau,omega,theta,m,R,L,I,xcorrido;
public:
  void Inicie(double theta0,double omega0, double m0,double R0,double L0,double x0corrido);
  void BorreFuerza(void);
  void AgregueFuerza(double tau0);
  void Mueva_theta(double dt,double Constante);
  void Mueva_omega(double dt,double Constante);
  void Dibujese(void);
  //  double Getx(void){return r.x();};
  //double GetFx(void){return r.y();};
  friend class Colisionador;
};

void Cuerpo::Inicie(double theta0,double omega0, double m0,double R0,double L0,double x0corrido){
  omega=omega0;
  theta=theta0;
  m=m0;
  R=R0;
  L=L0;
  xcorrido=x0corrido;
  I=m*L*L;
  
}

void Cuerpo::BorreFuerza(void){
  tau=0;
}

void Cuerpo::AgregueFuerza(double tau0){
  tau+=tau0;
}

void Cuerpo::Mueva_theta(double dt,double Constante){
theta+=omega*(Constante*dt);
}
void Cuerpo::Mueva_omega(double dt,double Constante){
omega+=tau*(Constante*dt)/I;
}

void Cuerpo::Dibujese(void){
  std::cout<<", "<<xcorrido+L*std::sin(theta)<<"+"<<R<<"*cos(t),"<<-L*std::cos(theta)<<"+"<<R<<"*sin(t)";
}
//------------------Funciones Globales---------

void InicieAnimacion(void){
  // std::cout<<"set terminal gif animate"<<std::endl;
  //std::cout<<"set output 'planeta.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set xrange [-12:22]"<<std::endl;
  std::cout<<"set yrange [-12:0]"<<std::endl;
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
 void  CalculeTodasLasFuerzas(Cuerpo* Pendulo);
 void  CalculeLaFuerzaEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2);
  
};

void Colisionador:: CalculeTodasLasFuerzas(Cuerpo* Pendulo){
  //Borrar todas las fuerzas
  for(int ii = 0;ii<N;ii++) Pendulo[ii].BorreFuerza();
  //Agregar fuerzas externas
  for(int ii = 0;ii<N;ii++) Pendulo[ii].AgregueFuerza(-Pendulo[ii].m*g*Pendulo[ii].L*std::sin(Pendulo[ii].theta));
  //calcular fuerza entre pares de pendulos
  for(int ii =0;ii<N;ii++){
    for(int jj =ii+1;jj<N;jj++){
      CalculeLaFuerzaEntre(Pendulo[ii],Pendulo[jj]);
    }
  }
  
}
void  Colisionador::CalculeLaFuerzaEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2){
  vector3D r[2];
  double L=Pendulo1.L;
  
  r[0].cargue(L*std::sin(Pendulo1.theta)+Pendulo1.xcorrido,-L*std::cos(Pendulo1.theta),0);
  r[1].cargue(L*std::sin(Pendulo2.theta)+Pendulo2.xcorrido,-L*std::cos(Pendulo2.theta),0);
    
  double r12=norma(r[0]-r[1]);
  double s = (Pendulo2.xcorrido-Pendulo1.xcorrido)-r12;
  double F;
  if(s>0){
    F=K*std::pow(s,1.5);
    Pendulo1.AgregueFuerza(-F*Pendulo1.L);
    Pendulo2.AgregueFuerza(F*Pendulo1.L);
  }
}



int main(void){
  
  double t, dt=0.001;
  double tdibujo;int Ndibujos;
  Cuerpo Pendulo[N]; 
  Colisionador Newton;
  
  double m0=1,L0=10;

  double R0=1;

  double x0corrido=2*R0;

  double theta0=-15*M_PI/180;

  double T= 2*M_PI*std::sqrt(L0/g), tmax=10*T;

  InicieAnimacion(); Ndibujos=2000;
  //---------------(theta0, omega0, m0, R0, L0, x0corrido);
  Pendulo[0].Inicie(theta0, 0     , m0, R0, L0, 0);
  for(int ii =1; ii<N;ii++){
    Pendulo[ii].Inicie(0, 0     , m0, R0, L0, x0corrido*ii);
  }
  
  Newton.CalculeTodasLasFuerzas(Pendulo);
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    if(tdibujo<tmax/Ndibujos){
      InicieCuadro();
      for(int ii = 0;ii<N;ii++)Pendulo[ii].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    

    
    //std::cout<<Pendulo[0].Getx()<<"   "<<Pendulo[0].Gety()<<std::endl;
    for(int ii = 0;ii<N;ii++)Pendulo[ii].Mueva_theta(dt,ZETA);
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(int ii = 0;ii<N;ii++)Pendulo[ii].Mueva_omega(dt,(1-2*LAMBDA)/2);
    for(int ii = 0;ii<N;ii++)Pendulo[ii].Mueva_theta(dt,CHI);
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(int ii = 0;ii<N;ii++)Pendulo[ii].Mueva_omega(dt,LAMBDA);
    for(int ii = 0;ii<N;ii++)Pendulo[ii].Mueva_theta(dt,1-2*(CHI+ZETA));
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(int ii = 0;ii<N;ii++)Pendulo[ii].Mueva_omega(dt,LAMBDA);		    
    for(int ii = 0;ii<N;ii++)Pendulo[ii].Mueva_theta(dt,CHI);
    Newton.CalculeTodasLasFuerzas(Pendulo);
    for(int ii = 0;ii<N;ii++)Pendulo[ii].Mueva_omega(dt,(1-2*LAMBDA)/2);
    for(int ii = 0;ii<N;ii++)Pendulo[ii].Mueva_theta(dt,ZETA);
  }
  
  return 0;
}
		    
  
