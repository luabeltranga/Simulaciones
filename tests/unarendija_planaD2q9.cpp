#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

const int Lx=1000;
const int Ly=600;

const int Q=9;
const double W0=4.0/9;



const double C=0.5;
const double TresC2=3*C*C;
const double AUX0=1-5*C*C/3.0;



const double tau=0.5;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

class LatticeBoltzmann{
private:
  double w[Q];
  int V[2][Q];  //V[alpha][i]  alpha=0 es x, alpha=1 es y
  double f[Lx][Ly][Q];
  double fnew[Lx][Ly][Q]; // f[ix][iy][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,bool UseNew,double sigma);
  double Jx(int ix,int iy);
  double Jy(int ix,int iy);
  double feq(int i,double rho0,double Jx0,double Jy0,int ix,int iy);
  void Inicie(double rho0,double Jx0,double Jy0);
  double GetSigma(int ix,int iy,int t);
  void Colisione(int t);
  void Adveccione(void);
  void Imprimase(char const * NombreArchivo, int t);
  void ImprimaseUnaLinea(char const * NombreArchivo, int t);
  void Max(int t);
  
};
LatticeBoltzmann::LatticeBoltzmann(void){
  w[0]=W0;
  w[1]=w[2]=w[3]=w[4]=1.0/9;
  w[5]=w[6]=w[7]=w[8]=1.0/36;

  V[0][0]=0;  
  V[1][0]=0;

  V[0][1]=1;    V[0][2]=0;    V[0][3]=-1;   V[0][4]=0;  
  V[1][1]=0;    V[1][2]=1;    V[1][3]=0;    V[1][4]=-1;

  V[0][5]=1;    V[0][6]=-1;    V[0][7]=-1;   V[0][8]=1;  
  V[1][5]=1;    V[1][6]=1;    V[1][7]=-1;    V[1][8]=-1;  
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
double LatticeBoltzmann::feq(int i,double rho0,double Jx0,double Jy0,int ix,int iy){
  if(i==0)
    return AUX0*rho0;
  else
    return w[i]*(TresC2*rho0+3*(V[0][i]*Jx0+V[1][i]*Jy0));
}

double LatticeBoltzmann::GetSigma(int ix,int iy,int t){
  double A=10,lambda=10,omega=2*M_PI*C/lambda;
  if(ix==0 && iy == 100)
    return A*sin(omega*t);
  else
    return 0;
}

void LatticeBoltzmann::Inicie(double rho0,double Jx0,double Jy0){
  int ix,iy,i;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++)
	f[ix][iy][i]=feq(i,rho0,Jx0,Jy0,ix,iy);
}
void LatticeBoltzmann::Colisione(int t){ //de f a fnew
  int ix,iy,i; double rho0,Jx0,Jy0; double sigma;
  for(ix=0;ix<Lx;ix++)
    for(iy=0;iy<Ly;iy++){ //Para cada celda
      sigma=GetSigma(ix,iy,t);
      if((ix>=50 && ix <55)){
	//primer hueco
	if(((iy>=0 && iy <=98) || (iy>=102 && iy<= 600))){
	  rho0=rho(ix,iy,false,sigma);  Jx0=0;  Jy0=0; //Calculo campos
	}
	else{
	  rho0=rho(ix,iy,false,sigma);  Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy); //Calculo campos
	}
      }
      else{
	rho0=rho(ix,iy,false,sigma);  Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy); //Calculo campos
      }
      for(i=0;i<Q;i++) //para cada dirección
	fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(i,rho0,Jx0,Jy0,ix,iy); //evoluciono
    }
}
void LatticeBoltzmann::Adveccione(void){ //de fnew a f
  int ix,iy,i,ixt,iyt;
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      for(i=0;i<Q;i++){
	ixt=(ix+V[0][i]+Lx)%Lx;
	iyt=(iy+V[1][i]+Ly)%Ly;
	f[ixt][iyt][i]=fnew[ix][iy][i];
      }
    }
  }
}
void LatticeBoltzmann::Imprimase(char const * NombreArchivo,int t){
  ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0; double sigma;
  for(int ix=0;ix<200;ix++){
    for(int iy=0;iy<Ly/3;iy++){
      sigma=GetSigma(ix,iy,t);
      rho0=rho(ix,iy,true,sigma);   Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy);
      //MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
      cout<<"splot [0:200] '-' "<<endl;
      cout<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    cout<<endl;
    //MiArchivo<<endl;
  }
  cout<<"e "<<endl;
  //MiArchivo.close();
}

void LatticeBoltzmann::ImprimaseUnaLinea(char const * NombreArchivo,int t){
  ofstream MiArchivo(NombreArchivo); double rho0,Jx0,Jy0; double sigma;
  int ix =150;
  for(int iy=0;iy<Ly/3;iy++){
      sigma=GetSigma(ix,iy,t);
      rho0=rho(ix,iy,true,sigma);   Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy);
      //MiArchivo<<ix<<" "<<iy<<" "<<rho0<<endl;
      cout<<"plot '-' w l"<<endl;
      cout<<iy<<" "<<rho0*rho0;
	
  }
  cout<<endl;
  //MiArchivo<<endl;
  cout<<"e "<<endl;
  //MiArchivo.close();
}

void LatticeBoltzmann::Max(int t){
  int ix =150, iy=100;
  double rho0,Jx0,Jy0,sigma;
  sigma=GetSigma(ix,iy,t);
  rho0=rho(ix,iy,true,sigma);   Jx0=Jx(ix,iy);  Jy0=Jy(ix,iy);
  cout<<t<<" "<<rho0*rho0;
  
  
  cout<<endl;
      
}


//---------------- Funciones Globales --------

int main(void){
  LatticeBoltzmann Ondas;
  int t,tmax=338;

  double rho0=0,Jx0=0,Jy0=0;

  //Inicie
  Ondas.Inicie(rho0,Jx0,Jy0);
  //cout<<"set pm3d map "<<endl;
  //cout<<"set size ratio 1 "<<endl;

  cout<<"set xrange [0:200] "<<endl;
  
  //Corra
  for(t=0;t<tmax;t++){
    Ondas.Colisione(t);
    Ondas.Adveccione();
    //Ondas.Imprimase("Espejo.dat",t);
    Ondas.ImprimaseUnaLinea("Espejo.dat",t);
    //Ondas.Max(t);
    //cerr<<t<<endl;
  }
  
  //Mostrar Resultado.

  
  return 0;
}
