#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"

const double K=1.0e4;
const double g=9.8;
const double Gamma=50;
const double Kcundall=10;
const double MU=0.4;
const double Lx=100, Ly=100;
const int Nx=5,Ny=5,N=Nx*Ny;


const double ZETA=0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double CHI=-0.06626458266981849;

class Cuerpo;
class Colisionador;

class Cuerpo{
private:
  vector3D r,V,F, omega, tau;double m,R,theta,I;
public:
  void Inicie(double x0,double y0,double z0,
	      double Vx0,double Vy0,double Vz0,
	      double theta0,double omega0,
	      double m0,double R0);
  void BorreFuerzayTorque(void);
  void AgregueFuerza(vector3D F0);
  void AgregueTorque(vector3D tau0);
  void Mueva_r(double dt,double Constante);
  void Mueva_V(double dt,double Constante);
  void Dibujese(void);
  double Getx(void){return r.x();};
  double Gety(void){return r.y();};
  friend class Colisionador;
};

void Cuerpo::Inicie(double x0,double y0,double z0,
		    double Vx0,double Vy0,double Vz0,double theta0,
		    double omega0,double m0,double R0){
  r.cargue(x0,y0,z0);
  V.cargue(Vx0,Vy0,Vz0);
  omega.cargue( 0, 0,omega0);
  theta=theta0;
  m=m0;
  R=R0;
  I=2.0/5.0*m*R*R;
}

void Cuerpo::BorreFuerzayTorque(void){
  F.cargue(0,0,0); tau.cargue(0,0,0);
}

void Cuerpo::AgregueTorque(vector3D tau0){
  tau+=tau0;
}
void Cuerpo::AgregueFuerza(vector3D F0){
  F+=F0;
}

void Cuerpo::Mueva_r(double dt,double Constante){
  r+=V*(Constante*dt); theta+=omega.z()*Constante*dt;
}

void Cuerpo::Mueva_V(double dt,double Constante){
  V+=F*(Constante*dt/m); omega+=tau*(Constante*dt/I);
}

void Cuerpo::Dibujese(void){
  std::cout<<", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t) , "
  <<r.x()<<"+"<<R*std::cos(theta)/7.0<<"*t,"<<r.y()<<"+"<<R*std::sin(theta)/7.0<<"*t";
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
  vector3D ele[N+4][N+4];bool EstabaEnColision[N+4][N+4];
  
public:
  void Inicie(void);
  void CalculeTodasLasFuerzas(Cuerpo* Grano,double dt);
  void CalculeLaFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2,
			    vector3D & ele, bool & EstabaEnColision, double dt);
};

void Colisionador:: Inicie(void){
  for(int ii=0;ii<N+4;ii++){
    for(int jj=ii+1;jj<N+4;jj++){
      ele[ii][jj].cargue(0,0,0); EstabaEnColision[ii][jj]=false;
    }
  }
}

void Colisionador:: CalculeTodasLasFuerzas(Cuerpo* Grano,double dt){
  vector3D g_vector; g_vector.cargue(0,-g,0);
  //Borrar todas las fuerzas y torques
  for(int ii = 0;ii<N+4;ii++) Grano[ii].BorreFuerzayTorque();

  //Fuerza de gravedad
  for(int ii = 0;ii<N;ii++) Grano[ii].AgregueFuerza(Grano[ii].m*g_vector);
  
  //Fuerza entre pares de bolas
  for(int ii =0;ii<N;ii++){
    for(int jj =ii+1;jj<N+4;jj++){
      CalculeLaFuerzaEntre(Grano[ii],Grano[jj],ele[ii][jj],EstabaEnColision[ii][jj],dt);
    }
  }
}
void  Colisionador::CalculeLaFuerzaEntre(Cuerpo & Grano1,Cuerpo & Grano2,vector3D & ele, bool & EstabaEnColision, double dt){
  vector3D r21, n, Vc, Vcn, Vct, t, Fn,Ft,F2 ; double R1, R2 ,d21 ,s ,m1 ,m2 ,m12, componenteVcn, normaVct, componenteFn,Ftmax,normaFt;
  double ERFF=1.0e-8;
  r21=Grano2.r-Grano1.r; d21=norma(r21);
  s=(Grano1.R+Grano2.R)- d21;
  
  if(s>0){
    //si se chocan
    //Geometria y dinamica del contacto
    m1=Grano1.m; m2=Grano2.m; m12=(m1*m2)/(m1+m2);
    R1=Grano1.R; R2=Grano2.R;
    n=r21/d21;
    //calcular velocidad de contacto y el vector tangente
    Vc=(Grano2.V-Grano1.V)-(Grano2.omega^n)*R2-(Grano1.omega^n)*R1;
    componenteVcn=Vc*n; Vcn=n*componenteVcn;Vct=Vc-Vcn; normaVct=norma(Vct);
    if(normaVct<ERFF) t.cargue(0,0,0); else t=Vct/normaVct;
	
    //Fuerza de hertz
    componenteFn=K*std::pow(s,1.5);
    //Disipacion plastica
    componenteFn-=m12*std::sqrt(s)*Gamma*componenteVcn;if(componenteFn<0) componenteFn=0;
    Fn=n*componenteFn;

    //fuerzas tangenciales
    //estatica
    ele+=(Vct*dt);
    Ft=ele*(-Kcundall);
    //cinetica
    Ftmax=MU*componenteFn;normaFt=norma(Ft);
    if(normaFt>Ftmax)Ft=ele*(-Ftmax/norma(ele));
    
      
    
    F2=Fn+Ft;
    Grano1.AgregueFuerza(F2*(-1));     Grano1.AgregueTorque((n*R1)^(Ft*(-1)));
    Grano2.AgregueFuerza(F2);          Grano2.AgregueTorque((n*(-R2))^Ft);
    
    EstabaEnColision=true;
  }
  else if(EstabaEnColision==true){
    ele.cargue(0,0,0);EstabaEnColision=false;
  }
  
}




int main(void){
  double t, dt=1.0e-3;
  double tdibujo;int Ndibujos;
  Cuerpo Grano[N+4]; 
  Colisionador Newton;
  Crandom ran64(1);double theta;
    
  double m0=1, R0=6,V=10,omega0=10;
  double Rpared=10000,Mpared=1000*m0;
  double dx=Lx/(Nx+1);double dy=Ly/(Ny+1);
  double T=Lx/V, tmax=3*T;

  InicieAnimacion(); Ndibujos=2000;
  //PAREDES
  //---------------(x0       , y0      ,z0 , Vx0, Vy0,Vz0 , m0    , R0);
  //Pared arriba
  Grano[N  ].Inicie(Lx/2     ,Ly+Rpared,  0, 0  , 0,0,0,0, Mpared, Rpared);
  //Pared abajo
  Grano[N+1].Inicie(Lx/2     ,-Rpared  ,  0, 0  , 0,0,0,0, Mpared, Rpared);
  //Pared derecha
  Grano[N+2].Inicie(Lx+Rpared, Ly/2    ,  0, 0  , 0,0,0,0, Mpared, Rpared);
  //Pared izquerda
  Grano[N+3].Inicie(-Rpared  , Ly/2    ,  0, 0  , 0,0,0,0, Mpared, Rpared);

  //GRANOS

  for(int ii=0;ii<Nx;ii++){
    for(int jj=0;jj<Ny;jj++){
      theta=2*M_PI*ran64.r();
      //--------------------(x0       , y0       ,z0, Vx0              , Vy0             ,Vz0, theta0, omega0, m0 , R0);
      Grano[ii+Nx*jj].Inicie((ii+1)*dx, (jj+1)*dy, 0, V*std::cos(theta),V*std::sin(theta),0  ,0      ,omega0 , m0 , R0);
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
		    
  
