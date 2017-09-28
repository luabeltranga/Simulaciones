#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"

const int L=4;
const int N=(2*L+1)*(2*L+1);
const double a = 50;


const double ZETA=0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double CHI=-0.06626458266981849;

class Nodos;
class Colisionador;

class Nodos{
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

void Nodos::Inicie(double x0,double y0,double z0,
		    double Vx0,double Vy0,double Vz0,double m0,double R0){
  r.cargue(x0,y0,z0);
  V.cargue(Vx0,Vy0,Vz0);
  m=m0;
  R=R0;
}

void Nodos::BorreFuerza(void){
  F.cargue(0,0,0); 
}

void Nodos::AgregueFuerza(vector3D F0){
  F+=F0;
}

void Nodos::Mueva_r(double dt,double Constante){
  r+=V*(Constante*dt); 
}

void Nodos::Mueva_V(double dt,double Constante){
  V+=F*(Constante*dt/m); 
}

void Nodos::Dibujese(void){
  std::cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}
//------------------Funciones Globales---------

void InicieAnimacion(void){
  //std::cout<<"set terminal gif animate"<<std::endl;
  //std::cout<<"set output 'planeta.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"unset border"<<std::endl;
  std::cout<<"unset xtics"<<std::endl;
  std::cout<<"unset ytics"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  std::cout<<"set xrange [-10:360]"<<std::endl;
  std::cout<<"set yrange [-10:670]"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange[0:15]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(void){
  std::cout<<"plot 0,0 ";
  std::cout<<" , "<<360/15<<"*t,0";
  std::cout<<" , "<<360/15<<"*t,660";
  std::cout<<" , 0,"<<660/15<<"*t";
  std::cout<<" , 360,"<<660/15<<"*t";
  
}
void TermineCuadro(void){
  std::cout<<std::endl;
}
//-------clase colisionador-------

class Colisionador{
private:
  
public:
  void CalculeTodasLasFuerzas(Nodos* Nodo,double dt);
  void CalculeLaFuerzaEntre(Nodos & Nodo1,Nodos & Nodo2);
			    
};


void Colisionador:: CalculeTodasLasFuerzas(Nodos* Nodo,double dt){
  
  //Borrar todas las fuerzas y torques
  for(int ii = 0;ii<N;ii++) Nodo[ii].BorreFuerza();
  /*
  //Fuerza entre pares de bolas
  for(int ii =0;ii<N;ii++){
    for(int jj =ii+1;jj<N+4;jj++){
      CalculeLaFuerzaEntre(Nodo[ii],Nodo[jj]);
    }
  }
  */
}
void  Colisionador::CalculeLaFuerzaEntre(Nodos & Nodo1,Nodos & Nodo2){


}




int main(void){
  double t, dt=1.0e-2;
  double tdibujo;int Ndibujos;
  double xran=0;
  double yran=0;
  
  Nodos Nodo[N]; 
  Colisionador Newton;
  Crandom ran64(1);
    
  double m0=1, R0=5;
  double tmax = 1000*dt;
  
  InicieAnimacion();
  Ndibujos=2000;
  //Malla base triangular regular
  
  for(int ii=0;ii<2*L+1;ii++){  

    for(int jj=0;jj<L;jj++){
      xran=(jj+1+0.5)*a;
      yran=ii*std::sqrt(3)*a;
      
      //------------------------(x0      , y0      ,z0, Vx0, Vy0, Vz0, m0 , R0);
      Nodo[ii*(2*L+1)+jj].Inicie(xran, yran, 0,   0,   0,  0 , m0 , R0);
    }
    for(int jj=L;jj<2*L+1;jj++){
      xran=(2*L-jj+1)*a;
      yran=(ii+1)*std::sqrt(3)*a*0.5;
     
      //----------------------------(x0  , y0  ,z0, Vx0, Vy0, Vz0, m0 , R0);
      Nodo[ii*(2*L+1)+jj].Inicie(xran, yran, 0,   0,   0,  0 , m0 , R0);
    }
  }

  
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){

    
    if(tdibujo>tmax/Ndibujos){
      InicieCuadro();
      for(int ii = 0;ii<N;ii++)Nodo[ii].Dibujese();
      TermineCuadro();
      tdibujo=0;
    }
    
    
    
    //std::cout<<Nodo[8].Getx()<<"   "<<Nodo[8].Gety()<<std::endl;
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_r(dt,ZETA);
    Newton.CalculeTodasLasFuerzas(Nodo,dt);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_V(dt,(1-2*LAMBDA)/2);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_r(dt,CHI);
    Newton.CalculeTodasLasFuerzas(Nodo,dt);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_V(dt,LAMBDA);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_r(dt,1-2*(CHI+ZETA));
    Newton.CalculeTodasLasFuerzas(Nodo,dt);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_V(dt,LAMBDA);		    
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_r(dt,CHI);
    Newton.CalculeTodasLasFuerzas(Nodo,dt);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_V(dt,(1-2*LAMBDA)/2);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_r(dt,ZETA);
  }
  
  return 0;
}
		    
  
