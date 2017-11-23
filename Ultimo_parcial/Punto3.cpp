#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=2000;
const int Ly=1;

const int Q=5;
const double W0=1.0/3;
const double D=0.6;

const double C=0.5; // C<0.707 celdas/click
const double TresC2=3*C*C;
const double AUX0=1-TresC2*(1-W0);

const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

class LatticeBoltzmann{
private:
  double w[Q];
  int V[2][Q];  //V[alpha][i]  alpha=0 es x, alpha=1 es y
  double f[Lx][Ly][Q],fnew[Lx][Ly][Q]; // f[ix][iy][i]
  double max[Lx];
  double min[Lx];
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,bool UseNew,double sigma);
  double Jx(int ix,int iy);
  double Jy(int ix,int iy);
  double feq(int i,double rho0,double Jx0,double Jy0);
  void Inicie(double rho0,double Jx0,double Jy0);
  double GetSigma(int ix,int iy,int t);
  void Colisione(int t);
  void Adveccione(void);
  void Envolventes(int t);
  double SWR(void);
  void Imprimase(char const * NombreArchivo,int t);
};
LatticeBoltzmann::LatticeBoltzmann(void){
  w[0]=W0;
  w[1]=w[2]=w[3]=w[4]=1.0/6;

  V[0][0]=0;  
  V[1][0]=0;

  V[0][1]=1;    V[0][2]=0;    V[0][3]=-1;   V[0][4]=0;  
  V[1][1]=0;    V[1][2]=1;    V[1][3]=0;    V[1][4]=-1;  
  for(int ix = 0;ix<Lx;ix++){
    max[ix]=-10.0;
    min[ix]=10.0;
  }
}
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew,double sigma){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
      suma+=fnew[ix][iy][i];
    else
      suma+=f[ix][iy][i];
  return suma+0.5*sigma;
}
double LatticeBoltzmann::Jx(int ix,int iy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    suma+=V[0][i]*f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::Jy(int ix,int iy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    suma+=V[1][i]*f[ix][iy][i];
  return suma;
}
double LatticeBoltzmann::feq(int i,double rho0,double Jx0,double Jy0){
  if(i==0)
    return AUX0*rho0;
  else
    return w[i]*(TresC2*rho0+3*(V[0][i]*Jx0+V[1][i]*Jy0));
}
double LatticeBoltzmann::GetSigma(int ix,int iy,int t){
  double A=1,lambda=1000,omega=2*M_PI*C/lambda;
  if(ix==0 )
    return A*sin(omega*t);
  else
    return 0;
}

void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0){
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++)
	f[ix][iy][i]=feq(i,rho0,Jx0,Jy0);
}
void LatticeBoltzmann::Colisione(int t){ //de f a fnew
  int ix,iy,i; double rho0,Jx0,Jy0; double sigma;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){ //Para cada celda
      sigma=GetSigma(ix,iy,t);
      rho0=rho(ix,iy,false,sigma);  Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy); 
      for(i=0;i<Q;i++){ //para cada dirección
	if(ix==Lx-1){
	  fnew[ix][iy][1]=D*f[ix][iy][3];
	  fnew[ix][iy][2]=D*f[ix][iy][4];
	  fnew[ix][iy][3]=D*f[ix][iy][1];
	  fnew[ix][iy][4]=D*f[ix][iy][2];
	}
	fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(i,rho0,Jx0,Jy0); 
      }
    }
}
void LatticeBoltzmann::Adveccione(void){ //de fnew a f
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++)
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i];
}
void LatticeBoltzmann::Imprimase(char const * NombreArchivo,int t){
  ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0; double sigma;
  for(int ix=1;ix<Lx-2;ix++){
    for(int iy=0;iy<Ly;iy++){
      MiArchivo<<ix<<" "<<max[ix]<<"   "<<min[ix];
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}

void LatticeBoltzmann::Envolventes(int t){
  double sigma,rho0;
  for(int ix = 0;ix<Lx;ix++){
    sigma=GetSigma(ix,0,t);
    rho0=rho(ix,0,true,sigma);   
    if(max[ix]<rho0)max[ix]=rho0;
    if(min[ix]>rho0)min[ix]=rho0;
  }
}

double LatticeBoltzmann::SWR(void){
  double test1 = -100;
  double test2 = 100;
  for(int ix = 0;ix<Lx;ix++){
    if(test1<max[ix])test1=max[ix];
    if(test2<min[ix])test2=min[ix];
  }
  return test1/test2;
}
//---------------- Funciones Globales --------

int main(void){
  LatticeBoltzmann Ondas;
  int t,tmax=84400;
  double SWR=M_PI, CReflexion=0;

  double rho0=0,Jx0=0,Jy0=0;

  //Inicie
  Ondas.Inicie(rho0,Jx0,Jy0);
  //Corra
  for(t=0;t<tmax;t++){
    Ondas.Colisione(t);
    Ondas.Adveccione();
    if(t>=80000)Ondas.Envolventes(t);
  }

  SWR=Ondas.SWR();
  CReflexion=pow((1-SWR)/(1+SWR),2);
  cout<<D<<" "<<CReflexion<<endl;
  return 0;
}
