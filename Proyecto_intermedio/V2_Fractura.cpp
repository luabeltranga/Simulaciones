//para correr en consola ./a.out "#depasos de tiempo" | gnuplot

#include <iostream>
#include <cmath>
#include "Vector.h"
#include "Random64.h"

const int L=20;
const int H=10;
const int N=(2*L+1)*H;
const double a = 50;
const double m0=0.86, R0=5;
const double EsobreA=10;
const double C = 3.0;


const double ZETA=0.1786178958448091;
const double LAMBDA=-0.2123418310626054;
const double CHI=-0.06626458266981849;

class Nodos;
class Resortes;

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
  friend class Resortes;
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
  std::cout<<"lc rgb \"black\", "<<r.x()<<"+"<<R<<"*cos(t),"<<r.y()<<"+"<<R<<"*sin(t)";
}

//-------clase Resortes-------

class Resortes{
private:
  bool   *Conectado =nullptr;
  double *Longitud_Natural=nullptr;
  double *ConstanteK=nullptr;
  double *Disorder=nullptr;
  double alpha = 0;
  bool Iniciar_Relajacion = true;
public:
  Resortes( Nodos * Nodo);
  ~Resortes(void);
  void Secar(void);
  void DibujeResortes(Nodos* Nodo);
  void CalculeTodasLasFuerzas(Nodos* Nodo,double dt);
  void CalculeLaFuerzaEntre(Nodos & Nodo1,Nodos & Nodo2,int ii,int jj);
			    
};

Resortes::Resortes(Nodos * Nodo){
  Longitud_Natural=new double [N*N];
  Conectado =new bool [N*N];
  ConstanteK= new double[N*N];
  Disorder= new double[N*N];
  double test =0;
  Crandom ran64(2309);
  for(int ii =0;ii<N;ii++){
    for(int jj = ii+1; jj<N;jj++){
      test=norma(Nodo[ii].r-Nodo[jj].r);
      if(test<1.2*a){
	Conectado[ii*N+jj]=true;
	Longitud_Natural[ii*N+jj]=test;
	ConstanteK[ii*N+jj]=test*EsobreA;
	Disorder[ii*N+jj]=std::fmod(ran64.r(),0.1);    
      }
      else{
	Conectado[ii*N+jj]=false;
	Longitud_Natural[ii*N+jj]=0.0;
	ConstanteK[ii*N+jj]=0.0;
	Disorder[ii*N+jj]=0.0;
      }
    }
  }

  //abajo
  for(int ii =0;ii<L;ii++){
    for(int jj =0;jj<L;jj++){
      Disorder[ii*N+jj]=10.0;
    }
  }

  //Arriba
  for(int ii =N-1;ii>N-L-2;ii--){
    for(int jj =N-1;jj>N-L-2;jj--){
      Disorder[ii*N+jj]=10.0;
    }
  }

  //lado derecho
  for(int ii =L-1;ii<N;){
    Disorder[ii*N+ii+1]=10.0;
    Disorder[(ii+1)*N+ii+2*L+1]=10.0;
    ii+=2*L+1;
  }
  //lado izquierdo
  for(int ii =0;ii<N;){
    Disorder[ii*N+ii+2*L]=10.0;
    Disorder[(ii+2*L)*N+ii+2*L+1]=10.0;
    ii+=2*L+1;
  }
 
}
Resortes::~Resortes(void){
  delete [] Conectado;
  delete [] Longitud_Natural;
  delete [] ConstanteK;
  delete [] Disorder;
}


void Resortes::Secar(void){
  int iitest=0, jjtest=0;
  double min=1.0;
  if(alpha<0.3){for(int ii =0;ii<N;ii++){
      for(int jj = 0; jj<N;jj++){
	if( Disorder[ii*N+jj]<min && (Disorder[ii*N+jj]!=0.0 && Disorder[ii*N+jj] <10.0)){
	iitest=ii;
	jjtest=jj;
	min=Disorder[ii*N+jj];
	}
      }
    }
  }
  for(int ii =0;ii<N;ii++){
    for(int jj = 0; jj<N;jj++){
      Longitud_Natural[ii*N+jj]=(1-min)*Longitud_Natural[ii*N+jj];
    }
  }
  
  Disorder[iitest*N+jjtest]=0.0;
  Conectado[iitest*N+jjtest]=false;
  Iniciar_Relajacion=true;
}
void Resortes::DibujeResortes(Nodos* Nodo){
  for(int ii =0;ii<N;ii++){
    for(int jj = ii+1; jj<N;jj++){
      
      if(Conectado[ii*N+jj]){
	std::cout<<std::showpos<<"lc rgb \"black\", "<<Nodo[ii].Getx()<<""<<(Nodo[jj].Getx()-Nodo[ii].Getx())/15.0<<"*t,"<<Nodo[ii].Gety()<<""<<(Nodo[jj].Gety()-Nodo[ii].Gety())/15<<"*t ";
      }
    }
  }
}
void Resortes:: CalculeTodasLasFuerzas(Nodos* Nodo,double dt){
  if(Iniciar_Relajacion){
  //Borrar todas las fuerzas 
  for(int ii = 0;ii<N;ii++) Nodo[ii].BorreFuerza();
  
  //Fuerza entre pares de bolas
  for(int ii = 0;ii<N;ii++){
    for(int jj =0;jj<N;jj++){
      if(Conectado[ii*N+jj])  CalculeLaFuerzaEntre(Nodo[ii],Nodo[jj],ii,jj);
    }
  }
  
 //Abajo
  for(int ii =0;ii<L;ii++){
    Nodo[ii].BorreFuerza();
  }

  //Arriba
  for(int ii =N-1;ii>N-L-2;ii--){
    Nodo[ii].BorreFuerza();
  }
  //lado derecho
  int k=0;
  for(int ii =1;ii<2*H;ii+=2){
    Nodo[ii*L+k].BorreFuerza(),Nodo[(ii+2)*L+k].BorreFuerza();
    k++;
  }

  //lado izquierdo
  k=0;
  for(int ii =0;ii<2*H;ii+=2){
    Nodo[ii*L+k].BorreFuerza(),Nodo[(ii+2)*L+k].BorreFuerza();
    k++;
  }
  
  }
  
}
void  Resortes::CalculeLaFuerzaEntre(Nodos & Nodo1,Nodos & Nodo2, int ii,int jj){
  vector3D r12,r_unitario,F12;
  double norma12 =0;
  double alpha =0;
  r12=Nodo2.r-Nodo1.r;
  norma12 = norma(r12);
  if(norma12 > 1.0e-9){
  r_unitario=r12*(1/norma12);
  F12=ConstanteK[ii*N+jj]*(norma12-a)*r_unitario;
  Nodo1.AgregueFuerza(F12-C*Nodo1.V); Nodo2.AgregueFuerza(F12*(-1)-C*Nodo2.V);
  alpha+=(Longitud_Natural[ii*N+jj]-norma12)/(Longitud_Natural[ii*N+jj]);
  if(alpha>Disorder[ii*N+jj]){
    Conectado[ii*N+jj]=false;
  }
  }
  else{
    F12.cargue(0,0,0);
    Nodo2.AgregueFuerza(F12); Nodo1.AgregueFuerza(F12*(-1));  
  }
}




//------------------Funciones Globales---------

void InicieAnimacion(void){
  //std::cout<<"set terminal gif animate"<<std::endl;
  //std::cout<<"set output 'planeta.gif'"<<std::endl;
  std::cout<<"unset key"<<std::endl;
  std::cout<<"unset border"<<std::endl;
  std::cout<<"unset xtics"<<std::endl;
  std::cout<<"unset ytics"<<std::endl;
  std::cout<<"set style line 5 lt rgb \"black\" lw 3 pt 6"<<std::endl;
  std::cout<<"set size ratio -1"<<std::endl;
  //std::cout<<"set xrange [-10:"<< (L+2)*a<<"]"<<std::endl;
  //std::cout<<"set yrange [-10:"<<1.5*H*a<<"]"<<std::endl;
  std::cout<<"set parametric"<<std::endl;
  std::cout<<"set trange[0:15]"<<std::endl;
  std::cout<<"set isosamples 12"<<std::endl;
}
void InicieCuadro(void){
  std::cout<<"plot 0,0 ";
    
}
void TermineCuadro(void){
  std::cout<<std::endl;
}

void Inicializa_malla(Nodos * Nodo){
 
  double xran=0,yran=0;
  double x0=0,y0=0;
  double R=0.1*a;
  Crandom ran64(1);

  double c=0,b=0,bold=0;

  for(int ii=0;ii<H;ii++){  
    for(int jj=0;jj<L;jj++){
      c=ran64.r();b=ran64.r();
      if(a>b){
	bold=b;
	b=c;
	c=bold;
      }
      x0=(jj+1+0.5)*a;  y0=(ii)*std::sqrt(3)*a;
      xran=(b*R*std::cos(2*M_PI*c/b))+x0;yran = (b*R*std::sin(2*M_PI*c/b))+y0;
      //------------------------(x0, y0,z0, Vx0, Vy0, Vz0, m0 , R0);
      Nodo[ii*(2*L+1)+jj].Inicie(xran, yran, 0,   0,   0,  0 , m0 , R0);
    }
    for(int jj=L;jj<2*L+1;jj++){
      c=ran64.r();b=ran64.r();
      if(a>b){
	bold=b;
	b=c;
	c=bold;
      }
      x0=(2*L-jj+1)*a;  y0=(ii+1-0.5)*std::sqrt(3)*a;
      xran=(b*R*std::cos(2*M_PI*a/c))+x0;yran = (b*R*std::sin(2*M_PI*c/b))+y0;
      //------------------------(x0, y0,z0, Vx0, Vy0, Vz0, m0 , R0);
      Nodo[ii*(2*L+1)+jj].Inicie(xran, yran, 0,   0,   0,  0 , m0 , R0);
    }
  } 
}


int main(int argc,char ** argv){

  int PasosDeTiempo = std::atoi(argv[1]);
  double t, dt=1.0e-2;
  double tdibujo;int Ndibujos;
  
  Nodos Nodo[N]; 
  Inicializa_malla(Nodo);
  Resortes Red (Nodo);

  int step=0;
    
  
  //  double tmax = 10*dt;
  double tmax = PasosDeTiempo*dt;
  
  InicieAnimacion();
  Ndibujos=500;
  //Malla base triangular regular
  
  
  for(t=tdibujo=0;t<tmax;t+=dt,tdibujo+=dt){
    step++;
    
    if(tdibujo>tmax/Ndibujos){
      InicieCuadro();
      for(int ii = 0;ii<N;ii++)Nodo[ii].Dibujese();
      Red.DibujeResortes(Nodo);
      TermineCuadro();
      tdibujo=0;
    }
    
    //if(step<10)Red.Secar();
    Red.Secar();
    
    //std::cout<<Nodo[8].Getx()<<"   "<<Nodo[8].Gety()<<std::endl;
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_r(dt,ZETA);
    Red.CalculeTodasLasFuerzas(Nodo,dt);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_V(dt,(1-2*LAMBDA)/2);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_r(dt,CHI);
    Red.CalculeTodasLasFuerzas(Nodo,dt);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_V(dt,LAMBDA);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_r(dt,1-2*(CHI+ZETA));
    Red.CalculeTodasLasFuerzas(Nodo,dt);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_V(dt,LAMBDA);		    
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_r(dt,CHI);
    Red.CalculeTodasLasFuerzas(Nodo,dt);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_V(dt,(1-2*LAMBDA)/2);
    for(int ii = 0;ii<N;ii++)Nodo[ii].Mueva_r(dt,ZETA);
    
  }
  
  return 0;
}
		    
  
