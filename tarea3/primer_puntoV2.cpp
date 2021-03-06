#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=400;
const int Ly=400;

const int Q=5;
const double W0=1.0/3;


const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

class LatticeBoltzmann{
private:
  double *w=nullptr;
  int *V=nullptr;
  double *f=nullptr;
  double *fnew=nullptr;
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  double rho(int ix,int iy,bool UseNew,double sigma);
  double Jx(int ix,int iy);
  double Jy(int ix,int iy);
  double feq(int i,double rho0,double Jx0,double Jy0,int ix,int iy);
  void Inicie(double rho0,double Jx0,double Jy0);
  double GetSigma(int ix,int iy,int t);
  double Ccelda(int ix,int iy);
  void Colisione(int t);
  void Adveccione(void);
  void Imprimase(char const * NombreArchivo,int t);
  void ImprimaUnaLinea(char const * NombreArchivo,int t);
  
};
LatticeBoltzmann::LatticeBoltzmann(void){
  w=new double[Q];
  V=new int[2*Q];  //V[alpha][i]  alpha=0 es x, alpha=1 es y
  f=new double[Lx*Ly*Q];
  fnew=new double[Lx*Ly*Q]; // f[ix*Ly*Q+iy*Q+i]
  w[0]=W0;
  w[1]=w[2]=w[3]=w[4]=1.0/6;

  V[0*2+0]=0;  
  V[1*2+0]=0;

  V[0*2+1]=1;    V[0*2+2]=0;    V[0*2+3]=-1;   V[0*2+4]=0;  
  V[1*2+1]=0;    V[1*2+2]=1;    V[1*2+3]=0;    V[1*2+4]=-1;  
}
LatticeBoltzmann::~LatticeBoltzmann(void){
  delete [] w;
  delete [] V;
  delete [] f;
  delete [] fnew;
}
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew,double sigma){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew)
      suma+=fnew[ix*Ly*Q+iy*Q+i];
    else
      suma+=f[ix*Ly*Q+iy*Q+i];
  return suma+0.5*sigma;
}
double LatticeBoltzmann::Ccelda(int ix,int iy){
  //double theta=60*M_PI/180;
  //return -0.125*std::tanh(sin(theta)*iy+cos(theta)*ix-100)+0.375;
  return 0.5;
}
double LatticeBoltzmann::Jx(int ix,int iy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    suma+=V[0*2+i]*f[ix*Ly*Q+iy*Q+i];
  return suma;
}
double LatticeBoltzmann::Jy(int ix,int iy){
  int i; double suma;
  for(suma=0,i=0;i<Q;i++)
    suma+=V[1*2+i]*f[ix*Ly*Q+iy*Q+i];
  return suma;
}
double LatticeBoltzmann::feq(int i,double rho0,double Jx0,double Jy0,int ix,int iy){
  double C = Ccelda(ix,iy);
  double TresC2=3*C*C;
  double AUX0=1-TresC2*(1-W0);
  
  if(i==0)
    return AUX0*rho0;
  else
    return w[i]*(TresC2*rho0+3*(V[0*2+i]*Jx0+V[1*2+i]*Jy0));
}

double LatticeBoltzmann::GetSigma(int ix,int iy,int t){
  double A=10,lambda=10,omega=2*M_PI*Ccelda(ix,iy)/lambda;
  if(ix==0)
    return A*sin(omega*t);
  else
    return 0;
}

void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0){
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++)
	f[ix*Ly*Q+iy*Q+i]=feq(i,rho0,Jx0,Jy0,ix,iy);
}
void LatticeBoltzmann::Colisione(int t){ //de f a fnew
  int ix,iy,i; double rho0,Jx0,Jy0; double sigma;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){ //Para cada celda
      sigma=GetSigma(ix,iy,t);
      rho0=rho(ix,iy,false,sigma);  Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy); //Calculo campos
      //	ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);     
	for(i=0;i<Q;i++) //para cada dirección
	  fnew[ix*Ly*Q+iy*Q+i]=UmUtau*f[ix*Ly*Q+iy*Q+i]+Utau*feq(i,rho0,Jx0,Jy0,ix,iy); //evoluciono
    }
}
void LatticeBoltzmann::Adveccione(void){ //de fnew a f
  int ix,iy,i,ixt,iyt;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(i=0;i<Q;i++){
	ixt=(ix+V[0*2+i]+Lx)%Lx;
	iyt=(iy+V[1*2+i]+Ly)%Ly;
	f[ixt*Ly*Q+iyt*Q+i]=fnew[ix*Ly*Q+iy*Q+i];
      }
    }
  }
}
void LatticeBoltzmann::Imprimase(char const * NombreArchivo,int t){
  ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0; double sigma;
  for(int ix=0;ix<Lx/2;ix++){
    for(int iy=0;iy<Ly;iy++){
      sigma=GetSigma(ix,iy,t);
      rho0=rho(ix,iy,true,sigma);   Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy);
      //      ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
      MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}
void LatticeBoltzmann::ImprimaUnaLinea(char const * NombreArchivo,int t){
  ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0; double sigma;
  int ix=Lx/2;
  for(int iy=0;iy<Ly;iy++){
    sigma=GetSigma(ix,iy,t);
    rho0=rho(ix,iy,true,sigma);   Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy);
    //    ImponerCampos(ix,iy,rho0,Jx0,Jy0,t);
    MiArchivo<<iy<<" "<<rho0<<endl;
  }
  MiArchivo.close();
}
//---------------- Funciones Globales --------

int main(void){
  LatticeBoltzmann Ondas;
  int t,tmax=600;

  double rho0=0,Jx0=0,Jy0=0;

  //Inicie
  Ondas.Inicie(rho0,Jx0,Jy0);
  //Corra
  for(t=0;t<tmax;t++){
    Ondas.Colisione(t);
    Ondas.Adveccione();
  }
  
  //Mostrar Resultado.
  Ondas.Imprimase("Ondas.dat",t);
  Ondas.ImprimaUnaLinea("CorteCentralConFuente.dat",t);

  return 0;
}
